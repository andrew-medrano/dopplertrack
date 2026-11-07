"""Candidate generation and scoring: start from the flights, not from a grid.

Two events heard dt apart at slant distances d1 and d2 mean the two aircraft were within
d1 + d2 of each other at emission times dt apart (give or take the propagation delays). That
is a cheap geometric join over a day's trajectories, and it leaves a few tens of thousands of
flight pairs, each with a small lens of ground where both could have been heard. Ground
points are laid only in those lenses and tested with exact closest approaches, gates and a
silence filter. A box is a restriction on the join, not a different search; a single event
with a known time is the same lens stage on one flight; a known location is the lens stage on
one point.

Stages, each a function of the previous one's output:
  1. join          pair_join: (segment 1, segment 2, bucket range) for every pair the
                   geometry allows, on the store's 10 s light index; single_join for one event
                   at a known time.
  2. lens          lens_candidates: per pair, a 1 km lattice gated loosely on the knot chords,
                   then a 200 m lattice gated exactly on the stored curves, then the silence
                   filter.
  3. tables        every candidate gets a first score by lookup in straight-pass score
                   tables (dopplertrack.tables).
  4. exact         the best `exact_points` candidates of each of the best `exact_pairs` pairs
                   by table score are scored on their real trajectories (locked_score).
  5. refine        a 50 m lattice around the best point of each of the best `refine_pairs`
                   pairs, gated and silence-filtered the same way, scored exactly.
  6. extras        rescore_with_extra: further flyovers on the top pairs, locked to the same
                   clock.
"""
from __future__ import annotations

import os
import sys
import time
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .cpa import arc_bound, closest_approach
from .doppler import speed_of_sound
from .geo import enu_to_lla, lla_to_enu
from .score import F0_SCALES, SHIFTS, EventWindow, locked_score, modelled_curve, single_score
from .store import TrackStore
from .tables import DELTAS, build_table

EXACT_PAIRS = 500        # pairs, by best table score, whose candidates are scored exactly
EXACT_POINTS = 20        # candidates per such pair, by table score
REFINE_PAIRS = 50        # pairs, by exact score, refined on a fine lattice
REFINE_HALF_M = 300.0
REFINE_STEP_M = 50.0
PRE_SLACK_M = 400.0      # distance slack of the fine lattice's pre-gate (chords, terrain height of the coarse cell)


def _log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def default_workers() -> int:
    """Process pool size: every core on Linux; 1 on macOS, where forking a process that has
    already used BLAS or Objective-C threads is unsafe. On macOS pass workers explicitly with
    OBJC_DISABLE_INITIALIZE_FORK_SAFETY=YES set in the environment."""
    return 1 if sys.platform == "darwin" else (os.cpu_count() or 1)


_G: dict = {}      # per-process state for the pools (set before forking)


def _pool_map(fn, items, workers: int, chunksize: int = 8):
    if workers <= 1:
        for it in items:
            yield fn(it)
        return
    from multiprocessing import get_context
    with get_context("fork").Pool(workers) as pool:
        yield from pool.imap_unordered(fn, items, chunksize=chunksize)


# ----------------------------------------------------------------------------------------
# the join
# ----------------------------------------------------------------------------------------

def lag_range(dt_heard: float, d1, d2, c: float, tol_s: float) -> tuple[float, float]:
    """Emission-time gap between the two CPAs allowed by the heard gap dt +- tol and the
    slant-distance gates: t_cpa2 - t_cpa1 = dt_heard + (d1 - d2) / c."""
    return dt_heard + (d1[0] - d2[1]) / c - tol_s, dt_heard + (d1[1] - d2[0]) / c + tol_s


def _empty_pairs():
    return pd.DataFrame({"seg1": np.zeros(0, np.int32), "seg2": np.zeros(0, np.int32), "b_min": np.zeros(0, np.int64),
                         "b_max": np.zeros(0, np.int64), "n": np.zeros(0, np.int64)})


def pair_join(st: TrackStore, dt_heard: float, v1, d1, v2, d2, tol_s: float = 6.0, c: float | None = None,
              v_slack: float = 0.1, a_mask: np.ndarray | None = None, progress=None) -> pd.DataFrame:
    """Flight pairs that could be heard as event 1 then event 2 by one observer: in some
    bucket k flight 1 (speed in v1) is within d1max + d2max of where flight 2 (speed in v2) is
    at the allowed lag. a_mask restricts the flight-1 samples (a bounding box). Returns one
    row per pair: seg1, seg2, b_min, b_max (buckets of the flight-1 side), n (matched
    buckets). Nothing here is exact; it only has to be complete."""
    c = speed_of_sound() if c is None else c
    W = st.bucket_s
    lo, hi = lag_range(dt_heard, d1, d2, c, tol_s)
    v1g = (v1[0] * (1 - v_slack), v1[1] * (1 + v_slack)); v2g = (v2[0] * (1 - v_slack), v2[1] * (1 + v_slack))
    R = d1[1] + d2[1] + (v1g[1] + v2g[1]) * W / 2
    sspd = np.asarray(st.sspd)
    A = (sspd >= v1g[0]) & (sspd <= v1g[1])
    if a_mask is not None:
        A &= a_mask
    B = (sspd >= v2g[0]) & (sspd <= v2g[1])
    buckets = np.unique(np.asarray(st.sb)[A])
    btrees: dict[int, tuple] = {}
    out1, out2, outb = [], [], []
    t0 = time.time()
    for n, k in enumerate(buckets):
        sl = st.bucket_slice(int(k)); a_idx = np.flatnonzero(A[sl]) + sl.start
        ta = cKDTree(np.asarray(st.spos[a_idx], float))
        tk = st.bucket_time(k)
        kb_lo, kb_hi = st.bucket(tk - W / 2 + lo), st.bucket(tk + W / 2 + hi)
        for kb in list(btrees):
            if kb < kb_lo:
                del btrees[kb]
        for kb in range(kb_lo, kb_hi + 1):
            if kb not in btrees:
                slb = st.bucket_slice(kb); b_idx = np.flatnonzero(B[slb]) + slb.start
                btrees[kb] = (b_idx, cKDTree(np.asarray(st.spos[b_idx], float))) if len(b_idx) else (b_idx, None)
            b_idx, tb = btrees[kb]
            if tb is None:
                continue
            m = ta.sparse_distance_matrix(tb, R, output_type="coo_matrix")
            if m.nnz == 0:
                continue
            s1 = st.ss[a_idx[m.row]]; s2 = st.ss[b_idx[m.col]]
            keep = s1 != s2
            out1.append(s1[keep]); out2.append(s2[keep]); outb.append(np.full(int(keep.sum()), k, np.int64))
        if progress and (n + 1) % 1000 == 0:
            progress(n + 1, len(buckets), time.time() - t0)
    if not out1:
        return _empty_pairs()
    df = pd.DataFrame({"seg1": np.concatenate(out1), "seg2": np.concatenate(out2), "b": np.concatenate(outb)}).drop_duplicates()
    g = df.groupby(["seg1", "seg2"], sort=False)["b"]
    return g.agg(b_min="min", b_max="max", n="size").reset_index()


def single_join(st: TrackStore, t_heard: float, v1, d1, tol_s: float = 30.0, c: float | None = None,
                v_slack: float = 0.1, a_mask: np.ndarray | None = None) -> pd.DataFrame:
    """One event heard at a known unix time (+- tol_s): every flight with a speed-gated sample
    in the buckets its CPA could fall in, as pair rows with seg2 = -1."""
    c = speed_of_sound() if c is None else c
    sspd = np.asarray(st.sspd); sb = np.asarray(st.sb)
    A = (sspd >= v1[0] * (1 - v_slack)) & (sspd <= v1[1] * (1 + v_slack))
    if a_mask is not None:
        A &= a_mask
    A &= (sb >= st.bucket(t_heard - d1[1] / c - tol_s)) & (sb <= st.bucket(t_heard + tol_s))
    if not A.any():
        return _empty_pairs()
    df = pd.DataFrame({"seg1": np.asarray(st.ss)[A], "b": sb[A]})
    g = df.groupby("seg1", sort=False)["b"].agg(b_min="min", b_max="max", n="size").reset_index()
    g.insert(1, "seg2", np.int32(-1))
    return g


# ----------------------------------------------------------------------------------------
# lens lattices, exact gates, silence filter
# ----------------------------------------------------------------------------------------

@dataclass
class LensParams:
    dt_heard: float                        # heard gap event 2 - event 1 (two events) or the event's heard unix time (single)
    v1: tuple; d1: tuple
    v2: tuple = (0.0, 0.0); d2: tuple = (0.0, 0.0)
    tol_s: float = 6.0
    c: float = 340.0
    spacing_m: float = 200.0
    coarse_m: float = 1000.0
    max_range_m: float = 12000.0
    silence_d_m: float | None = 6000.0
    silence_margin_s: float = 40.0
    pad_s: float = 30.0
    single: bool = False                   # one event at a known time: dt_heard is that time
    points: tuple | None = None            # (lat, lon, elev, enu (G,3)) to test instead of a lattice

    def fingerprint(self) -> dict:
        out = {}
        for k, v in self.__dict__.items():
            if k == "points":
                out[k] = None if v is None else [float(x) for x in np.asarray(v[0]).ravel()[:50]]
            else:
                out[k] = tuple(float(x) for x in v) if isinstance(v, (tuple, list)) else v
        return out

    @property
    def d1hi(self) -> float:
        return min(self.d1[1], self.max_range_m)

    @property
    def d2hi(self) -> float:
        return min(self.d2[1], self.max_range_m)


def _ground_lattice(e, n, lat_c, lon_c, elevation, lat0, lon0):
    """Local tangent-plane offsets (e, n) at (lat_c, lon_c) -> lat, lon, elev, global ENU."""
    lat, lon, _ = enu_to_lla(e, n, np.zeros_like(e), lat_c, lon_c)
    el = np.maximum(elevation(lat, lon), 0.0) if elevation is not None else np.zeros_like(lat)
    return lat, lon, el, lla_to_enu(lat, lon, el, lat0, lon0)


def _gate(G, ch1, ch2, p: LensParams, exact=True, slack_d=0.0, slack_t=0.0, v_slack=0.0, dist_only=False):
    """The exact gates on ground points G: CPA distance, speed and interior on flight 1 (and on
    flight 2), and the heard gap between the two CPAs (or the heard time of the one event).
    Flight 2's closest approach is computed only for the points flight 1 lets through.
    dist_only applies the distance gates alone: a pre-gate on chords must not judge time or
    speed, since on a turning track the chord step's CPA time and speed can be off by more
    than any small slack (its distance cannot: it is within the arc bound).
    Returns (ok, da, ha, va, db, hb, vb); flight-2 values are NaN where not computed."""
    da, ha, va, ia = closest_approach(G, ch1, p.c, exact)
    ok = (da >= p.d1[0] - slack_d) & (da <= p.d1hi + slack_d)
    if not dist_only:
        ok &= ia & (va >= p.v1[0] * (1 - v_slack)) & (va <= p.v1[1] * (1 + v_slack))
        if p.single:
            ok &= np.abs(ha - p.dt_heard) <= p.tol_s + slack_t
    db, hb, vb = (np.full(len(G), np.nan) for _ in range(3))
    if ch2 is not None and ok.any():
        i = np.flatnonzero(ok)
        d2, h2, v2, i2 = closest_approach(G[i], ch2, p.c, exact)
        ok2 = (d2 >= p.d2[0] - slack_d) & (d2 <= p.d2hi + slack_d)
        if not dist_only:
            ok2 &= (i2 & (v2 >= p.v2[0] * (1 - v_slack)) & (v2 <= p.v2[1] * (1 + v_slack))
                    & (np.abs(h2 - ha[i] - p.dt_heard) <= p.tol_s + slack_t))
        ok[i] = ok2; db[i], hb[i], vb[i] = d2, h2, v2
    return ok, da, ha, va, db, hb, vb


def _silent(st: TrackStore, p: LensParams, F, ha, hb, seg1: int, seg2: int) -> np.ndarray:
    """True where no other flight passes within silence_d_m of the point with its CPA heard
    between the two events (by more than silence_margin_s): a third flyover in between would
    have been heard. Points are taken in groups by their event-1 time (a pair whose flights
    stay together for an hour has candidates whose event windows are an hour apart), each
    group against the flights near it in its own window."""
    keep = np.ones(len(F), bool)
    if p.silence_d_m is None or p.single or not len(F):
        return keep
    span = 600.0
    grp = np.floor((ha - ha.min()) / span).astype(np.int64)
    for gi in np.unique(grp):
        i = np.flatnonzero(grp == gi)
        keep[i] = _silent_group(st, p, F[i], ha[i], hb[i], seg1, seg2)
    return keep


def _silent_group(st: TrackStore, p: LensParams, F, ha, hb, seg1: int, seg2: int) -> np.ndarray:
    keep = np.ones(len(F), bool)
    W, mg = st.bucket_s, p.silence_margin_s
    ta, tb = float(ha.min()) + mg, float(hb.max()) - mg
    if tb <= ta:
        return keep
    ctr = F.mean(0); rad = float(np.linalg.norm(F - ctr, axis=1).max())
    others: set[int] = set()
    for b in range(st.bucket(ta), st.bucket(tb) + 1):
        tr = st.tree(b)
        if tr is None:
            continue
        sl = st.bucket_slice(b)
        hit = tr.query_ball_point(ctr, rad + p.silence_d_m + 150.0 * W / 2)
        others.update(int(x) for x in st.ss[np.asarray(hit, int) + sl.start])
    others.discard(seg1); others.discard(seg2)
    if not others:
        return keep
    tf = cKDTree(F)
    for so in sorted(others):
        if not keep.any():
            break
        ch = st.chain(so, ta - 60.0, tb + 60.0)
        if ch is None:
            continue
        # only points within silence_d_m of the curve can be silenced: those within that plus
        # the longest half-chord plus the arc bound of some chord midpoint
        mid = 0.5 * (ch.p[1:] + ch.p[:-1])
        reach = p.silence_d_m + 0.5 * float(np.linalg.norm(ch.p[1:] - ch.p[:-1], axis=1).max()) + arc_bound(ch)
        hits = tf.query_ball_point(mid, reach)
        near = np.unique(np.concatenate([np.asarray(h, np.int64) for h in hits])) if len(hits) else np.zeros(0, np.int64)
        near = near[keep[near]]          # a point already silenced stays silenced
        if not len(near):
            continue
        do, ho, _, io = closest_approach(F[near], ch, p.c)
        keep[near] &= ~(io & (do <= p.silence_d_m) & (ho > ha[near] + mg) & (ho < hb[near] - mg))
    return keep


def _frame(lat, lon, el, F, seg1, seg2, ha, da, va, hb, db, vb) -> pd.DataFrame:
    return pd.DataFrame({"lat": lat, "lon": lon, "elev": el, "e": F[:, 0], "n": F[:, 1], "u": F[:, 2],
                         "seg1": np.int32(seg1), "seg2": np.int32(seg2), "t_heard1": ha, "d1": da.astype(np.float32),
                         "v1": va.astype(np.float32), "t_heard2": hb, "d2": db.astype(np.float32), "v2": vb.astype(np.float32)})


def _lens_one(row) -> pd.DataFrame | None:
    st: TrackStore = _G["store"]; p: LensParams = _G["params"]; elevation = _G.get("elevation")
    seg1, seg2, b_min, b_max = int(row[0]), int(row[1]), int(row[2]), int(row[3])
    t1a, t1b = st.bucket_start(b_min) - p.pad_s, st.bucket_start(b_max + 1) + p.pad_s
    ch1 = st.chain(seg1, t1a, t1b)
    if ch1 is None:
        return None
    ch2 = None
    if not p.single:
        lo, hi = lag_range(p.dt_heard, p.d1, p.d2, p.c, p.tol_s)
        ch2 = st.chain(seg2, t1a + lo - p.pad_s, t1b + hi + p.pad_s)
        if ch2 is None:
            return None
    if p.points is not None:
        lat, lon, el, G = (np.asarray(x, float) for x in p.points)
    else:
        # coarse lattice in the tangent plane at the flight-1 window's centre, within d1hi of its
        # ground track (the knot chain, densified so the distance test sees the whole path)
        lat_c, lon_c, _ = enu_to_lla(*ch1.p.mean(0), st.lat0, st.lon0)
        lat1, lon1, _ = enu_to_lla(ch1.p[:, 0], ch1.p[:, 1], ch1.p[:, 2], st.lat0, st.lon0)
        P1 = lla_to_enu(lat1, lon1, np.zeros_like(lat1), lat_c, lon_c)[:, :2]
        u = np.linspace(0, 1, 8, endpoint=False)[None, :, None]
        P1 = np.concatenate([(P1[:-1, None, :] + u * (P1[1:] - P1[:-1])[:, None, :]).reshape(-1, 2), P1[-1:]])
        cm, d1hi = p.coarse_m, p.d1hi
        e = np.arange(np.floor((P1[:, 0].min() - d1hi) / cm), np.ceil((P1[:, 0].max() + d1hi) / cm) + 1) * cm
        n = np.arange(np.floor((P1[:, 1].min() - d1hi) / cm), np.ceil((P1[:, 1].max() + d1hi) / cm) + 1) * cm
        E, N = np.meshgrid(e, n, indexing="ij"); E, N = E.ravel(), N.ravel()
        near = cKDTree(P1).query(np.stack([E, N], 1), distance_upper_bound=d1hi + cm)[0] <= d1hi + cm
        E, N = E[near], N[near]
        if len(E) == 0:
            return None
        _, _, _, G = _ground_lattice(E, N, lat_c, lon_c, None, st.lat0, st.lon0)
        sd = 0.75 * cm; stt = sd * (1.0 / max(p.v1[0] if p.single else min(p.v1[0], p.v2[0]), 10.0) + 1.0 / p.c)
        ok = _gate(G, ch1, ch2, p, slack_d=sd, slack_t=2 * stt, v_slack=0.1)[0]
        if not ok.any():
            return None
        # fine lattice around each surviving coarse cell; each fine point first takes its coarse
        # cell's terrain height and is gated with slack, and only the survivors look up their own
        Ec, Nc = E[ok], N[ok]
        elc = (np.maximum(elevation(*_ground_lattice(Ec, Nc, lat_c, lon_c, None, st.lat0, st.lon0)[:2]), 0.0)
               if elevation is not None else np.zeros(len(Ec)))
        m = int(np.ceil(cm / 2 / p.spacing_m)); offs = np.arange(-m, m + 1) * p.spacing_m
        OE, ON = np.meshgrid(offs, offs, indexing="ij")
        Ef = (Ec[:, None] + OE.ravel()[None, :]).ravel(); Nf = (Nc[:, None] + ON.ravel()[None, :]).ravel()
        Elf = np.repeat(elc, OE.size)
        key, first = np.unique(np.rint(np.stack([Ef, Nf], 1) / p.spacing_m).astype(np.int64), axis=0, return_index=True)
        Ef, Nf, Elf = key[:, 0] * p.spacing_m, key[:, 1] * p.spacing_m, Elf[first]
        la0, lo0, _ = enu_to_lla(Ef, Nf, np.zeros_like(Ef), lat_c, lon_c)
        G0 = lla_to_enu(la0, lo0, Elf, st.lat0, st.lon0)
        pre = _gate(G0, ch1, ch2, p, exact=False, slack_d=PRE_SLACK_M, dist_only=True)[0]
        if not pre.any():
            return None
        lat, lon, el, G = _ground_lattice(Ef[pre], Nf[pre], lat_c, lon_c, elevation, st.lat0, st.lon0)
    ok, da, ha, va, db, hb, vb = _gate(G, ch1, ch2, p)
    if not ok.any():
        return None
    F = G[ok]; ha, da, va, hb, db, vb = ha[ok], da[ok], va[ok], hb[ok], db[ok], vb[ok]
    lat, lon, el = lat[ok], lon[ok], el[ok]
    keep = _silent(st, p, F, ha, hb, seg1, seg2)
    if not keep.any():
        return None
    F, ha, da, va, hb, db, vb, lat, lon, el = (x[keep] for x in (F, ha, da, va, hb, db, vb, lat, lon, el))
    return _frame(lat, lon, el, F, seg1, seg2, ha, da, va, hb, db, vb)


CAND_COLUMNS = ["lat", "lon", "elev", "e", "n", "u", "seg1", "seg2", "t_heard1", "d1", "v1", "t_heard2", "d2", "v2"]


def empty_candidates() -> pd.DataFrame:
    return pd.DataFrame({k: np.zeros(0, np.int32 if k in ("seg1", "seg2") else float) for k in CAND_COLUMNS})


def lens_candidates(st: TrackStore, pairs: pd.DataFrame, params: LensParams, elevation=None,
                    workers: int | None = None, progress=None) -> pd.DataFrame:
    """Exact (point, flight pair) candidates for every pair from pair_join / single_join: a
    200 m lattice (params.spacing_m) wherever a 1 km (params.coarse_m) sweep passed loosened
    gates, with the gates and silence filter applied exactly; or just params.points."""
    workers = default_workers() if workers is None else workers
    _G.update(store=st, params=params, elevation=elevation)
    rows = pairs[["seg1", "seg2", "b_min", "b_max"]].to_numpy()
    out = []
    t0 = time.time()
    for k, df in enumerate(_pool_map(_lens_one, list(rows), workers, chunksize=64)):
        if df is not None:
            out.append(df)
        if progress and (k + 1) % 2000 == 0:
            progress(k + 1, len(rows), time.time() - t0)
    return pd.concat(out, ignore_index=True)[CAND_COLUMNS] if out else empty_candidates()


def save_candidates(path, cands: pd.DataFrame, fingerprint: dict):
    np.savez_compressed(path, fingerprint=np.array(repr(fingerprint)), **{k: cands[k].to_numpy() for k in CAND_COLUMNS})


def load_candidates(path, fingerprint: dict) -> pd.DataFrame | None:
    z = np.load(path, allow_pickle=False)
    if str(z["fingerprint"]) != repr(fingerprint):
        return None
    return pd.DataFrame({k: z[k] for k in CAND_COLUMNS})


# ----------------------------------------------------------------------------------------
# exact scoring
# ----------------------------------------------------------------------------------------

RESULT_COLUMNS = ["gidx", "lat", "lon", "elev", "seg1", "seg2", "icao1", "icao2", "t_heard1", "d1", "v1", "d2", "v2",
                  "score", "s1", "s2", "shift", "stage", "lat0", "lon0"]


def _score_chunk(args):
    st: TrackStore = _G["store"]; w1: EventWindow = _G["w1"]; w2 = _G["w2"]; c = _G["c"]
    i0, seg1, seg2, obs, th1 = args
    out = np.empty((len(seg1), 4))
    for k in range(len(seg1)):
        t0_map = w1.fit.t_cpa - float(th1[k])
        tr1 = st.track(int(seg1[k]), *w1.emission_span(t0_map))
        if seg2[k] < 0:
            s, sh = single_score(w1, tr1, obs[k], float(th1[k]), c=c)
            out[k] = (s, s, 0.0, sh)
        else:
            tr2 = st.track(int(seg2[k]), *w2.emission_span(t0_map))
            out[k] = locked_score(w1, w2, tr1, tr2, obs[k], float(th1[k]), c=c)
    return i0, out


def score_candidates(w1: EventWindow, w2: EventWindow | None, st: TrackStore, cands: pd.DataFrame, c: float,
                     workers: int | None = None, chunk: int = 100, progress=None) -> np.ndarray:
    """locked_score for every candidate row (the single-event sweep where seg2 < 0), on the
    stored trajectories sampled at 1 Hz. Returns (n, 4): score, s1, s2, shift."""
    workers = default_workers() if workers is None else workers
    w1.comb(); w2 is not None and w2.comb()                # built once, before any fork
    _G.update(store=st, w1=w1, w2=w2, c=c)
    obs = cands[["e", "n", "u"]].to_numpy(); s1 = cands.seg1.to_numpy(); s2 = cands.seg2.to_numpy(); th = cands.t_heard1.to_numpy()
    jobs = [(a, s1[a:a + chunk], s2[a:a + chunk], obs[a:a + chunk], th[a:a + chunk]) for a in range(0, len(cands), chunk)]
    sc = np.zeros((len(cands), 4))
    t0 = time.time(); done = 0
    for i0, out in _pool_map(_score_chunk, jobs, workers, chunksize=1):
        sc[i0:i0 + len(out)] = out; done += len(out)
        if progress:
            progress(done, len(cands), time.time() - t0)
    return sc


def _results(st: TrackStore, cands: pd.DataFrame, sc: np.ndarray, stage) -> pd.DataFrame:
    s1 = cands.seg1.to_numpy(); s2 = cands.seg2.to_numpy()
    res = pd.DataFrame({
        "gidx": np.arange(len(cands)), "lat": cands.lat.to_numpy(), "lon": cands.lon.to_numpy(), "elev": cands.elev.to_numpy(),
        "seg1": s1, "seg2": s2, "icao1": st.icao[s1], "icao2": np.where(s2 >= 0, st.icao[np.maximum(s2, 0)], ""),
        "t_heard1": cands.t_heard1.to_numpy(), "d1": cands.d1.to_numpy(), "v1": cands.v1.to_numpy(),
        "d2": cands.d2.to_numpy(), "v2": cands.v2.to_numpy(),
        "score": sc[:, 0], "s1": sc[:, 1], "s2": sc[:, 2], "shift": sc[:, 3], "stage": stage})
    res["lat0"] = st.lat0; res["lon0"] = st.lon0
    return res


def best_per_pair(res: pd.DataFrame) -> pd.DataFrame:
    """One row per flight pair (its best point), best pair first."""
    return res.sort_values("score", ascending=False).drop_duplicates(["seg1", "seg2"]).reset_index(drop=True)


# ----------------------------------------------------------------------------------------
# tables -> shortlist, and the 50 m refinement
# ----------------------------------------------------------------------------------------

def table_scores(cands: pd.DataFrame, t1, t2, dt_heard: float) -> np.ndarray:
    """First score of every candidate from the straight-pass tables: event 1 at its best shift,
    event 2 at the offset its geometry and that shift put it (as locked_score places it)."""
    a1 = t1.lookup(cands.d1.to_numpy(), cands.v1.to_numpy())
    k = np.argmax(a1, 1); s = a1[np.arange(len(k)), k].astype(float)
    if t2 is not None:
        delta = (cands.t_heard2.to_numpy() - cands.t_heard1.to_numpy()) - dt_heard + t1.shifts[k]
        s = s + t2.lookup_at(cands.d2.to_numpy(), cands.v2.to_numpy(), delta)
    return s


def shortlist(cands: pd.DataFrame, score: np.ndarray, top_pairs: int, top_points: int) -> np.ndarray:
    """Row indices of the best top_points candidates (by score) of each of the best top_pairs
    pairs (by their best candidate)."""
    df = pd.DataFrame({"seg1": cands.seg1.to_numpy(), "seg2": cands.seg2.to_numpy(), "s": score})
    best = df.groupby(["seg1", "seg2"], sort=False).s.max().sort_values(ascending=False).head(top_pairs)
    keep = df.set_index(["seg1", "seg2"]).index.isin(best.index)
    rank = df.groupby(["seg1", "seg2"], sort=False).s.rank(ascending=False, method="first").to_numpy()
    return np.flatnonzero(keep & (rank <= top_points))


def refine_candidates(st: TrackStore, res: pd.DataFrame, p: LensParams, elevation=None, top_pairs: int = REFINE_PAIRS,
                      half_m: float = REFINE_HALF_M, step_m: float = REFINE_STEP_M) -> pd.DataFrame:
    """A step_m lattice of half-width half_m around the best point of each of the best
    top_pairs pairs of an exactly scored result, gated and silence-filtered exactly like the
    lens (the centre point itself is left out: it is already scored)."""
    b = best_per_pair(res).head(top_pairs)
    off = np.arange(-half_m, half_m + step_m / 2, step_m)
    OE, ON = np.meshgrid(off, off, indexing="ij")
    centre = (np.abs(OE) < 1e-9) & (np.abs(ON) < 1e-9)
    OE, ON = OE[~centre], ON[~centre]
    out = []
    for row in b.itertuples(index=False):
        s1, s2 = int(row.seg1), int(row.seg2)
        lat, lon, el, G = _ground_lattice(OE, ON, row.lat, row.lon, elevation, st.lat0, st.lon0)
        te1 = float(row.t_heard1) - float(row.d1) / p.c
        ch1 = st.chain(s1, te1 - 150.0, te1 + 150.0)
        ch2 = None
        if s2 >= 0:
            te2 = te1 + p.dt_heard + (float(row.d1) - float(row.d2)) / p.c
            ch2 = st.chain(s2, te2 - 150.0, te2 + 150.0)
            if ch2 is None:
                continue
        if ch1 is None:
            continue
        ok, da, ha, va, db, hb, vb = _gate(G, ch1, ch2, p)
        if not ok.any():
            continue
        F = G[ok]; ha, da, va, hb, db, vb, lat_, lon_, el_ = (x[ok] for x in (ha, da, va, hb, db, vb, lat, lon, el))
        keep = _silent(st, p, F, ha, hb, s1, s2)
        if keep.any():
            out.append(_frame(lat_[keep], lon_[keep], el_[keep], F[keep], s1, s2, ha[keep], da[keep], va[keep],
                              hb[keep], db[keep], vb[keep]))
    return pd.concat(out, ignore_index=True)[CAND_COLUMNS] if out else empty_candidates()


# ----------------------------------------------------------------------------------------
# the whole search
# ----------------------------------------------------------------------------------------

def run_search(sp, ev1, ev2, st: TrackStore, harmonics=(1, 1), tol_s=None, v_tol=0.35, d_tol=(0.4, 2.5),
               silence_d_m=6000.0, max_range_m=12000.0, spacing_m=200.0, coarse_m=1000.0, elevation=None,
               c=None, score_mode="brightness", workers=None, cand_cache=None, bbox=None, points=None,
               t_heard_unix=None, exact_pairs: int = EXACT_PAIRS, exact_points: int = EXACT_POINTS,
               refine_pairs: int = REFINE_PAIRS, refine_half_m: float = REFINE_HALF_M,
               refine_step_m: float = REFINE_STEP_M, log=_log) -> pd.DataFrame:
    """Join -> lens candidates -> table scores -> exact scores of the shortlist -> 50 m
    refinement, on one track store. Returns the exactly scored rows (RESULT_COLUMNS, `stage`
    "exact" or "refine"), best first; .attrs["stats"] holds counts and stage times (s).

    ev2=None with t_heard_unix: one event at a known time (tol_s defaults to 30 s).
    bbox: restrict the join to flights in (lat_min, lat_max, lon_min, lon_max).
    points: (lat, lon, elev) arrays to test instead of lattices (a known location); every
    candidate is then scored exactly and nothing is refined.
    cand_cache: optional .npz path; the candidate table depends on the store, the events'
    gates and the geometry parameters but not on the spectrogram or scorer, so it is reused
    when those match.
    When no pair has more than exact_points candidates and there are at most exact_pairs
    pairs, the tables are skipped and every candidate is scored exactly."""
    c = speed_of_sound() if c is None else c
    workers = default_workers() if workers is None else workers
    single = ev2 is None
    if single and t_heard_unix is None:
        raise ValueError("one event needs its heard time (t_heard_unix)")
    tol_s = (30.0 if single else 6.0) if tol_s is None else tol_s
    g = lambda ev: ((ev.v * (1 - v_tol), ev.v * (1 + v_tol)), (ev.d * d_tol[0], min(ev.d * d_tol[1], max_range_m)))
    v1, d1 = g(ev1)
    v2, d2 = g(ev2) if not single else ((0.0, 0.0), (0.0, 0.0))
    pts = None
    if points is not None:
        la, lo, el = (np.atleast_1d(np.asarray(x, float)) for x in points)
        pts = (la, lo, el, lla_to_enu(la, lo, el, st.lat0, st.lon0))
    params = LensParams(dt_heard=float(t_heard_unix if single else ev2.t_cpa - ev1.t_cpa), v1=v1, d1=d1, v2=v2, d2=d2,
                        tol_s=tol_s, c=float(c), spacing_m=spacing_m, coarse_m=coarse_m, max_range_m=max_range_m,
                        silence_d_m=silence_d_m, single=single, points=pts)
    T: dict[str, float] = {}
    stats: dict = {"times": T}
    fp = {"store": st.fingerprint(), "params": params.fingerprint(), "dem": elevation is not None,
          "bbox": None if bbox is None else tuple(float(x) for x in bbox)}
    cands = None
    if cand_cache and os.path.exists(cand_cache):
        cands = load_candidates(cand_cache, fp)
        log(f"loaded {len(cands)} candidates from {cand_cache}" if cands is not None else f"{cand_cache} is for other inputs; recomputing")
    if cands is None:
        t0 = time.time()
        a_mask = st.sample_mask_in_bbox(bbox, pad_m=d1[1]) if bbox is not None else None
        if pts is not None:
            near = cKDTree(pts[3]).query(np.asarray(st.spos, float), distance_upper_bound=d1[1] + 1500.0)[0] <= d1[1] + 1500.0
            a_mask = near if a_mask is None else a_mask & near
        if single:
            pairs = single_join(st, params.dt_heard, v1, d1, tol_s=tol_s, c=c, a_mask=a_mask)
            log(f"join: {len(pairs)} flights near the event time pass the speed gate; gates v {v1[0]:.0f}-{v1[1]:.0f} m/s, d {d1[0]:.0f}-{d1[1]:.0f} m")
        else:
            pairs = pair_join(st, params.dt_heard, v1, d1, v2, d2, tol_s=tol_s, c=c, a_mask=a_mask,
                              progress=lambda k, n, s: log(f"  join: bucket {k}/{n}, {s:.0f} s"))
            log(f"join: {len(pairs)} flight pairs; gates v1 {v1[0]:.0f}-{v1[1]:.0f} m/s, d1 {d1[0]:.0f}-{d1[1]:.0f} m, "
                f"v2 {v2[0]:.0f}-{v2[1]:.0f} m/s, d2 {d2[0]:.0f}-{d2[1]:.0f} m")
        T["join"] = time.time() - t0
        stats["pairs_joined"] = len(pairs)
        t0 = time.time()
        cands = lens_candidates(st, pairs, params, elevation=elevation, workers=workers,
                                progress=lambda k, n, s: log(f"  lens: pair {k}/{n}, {s:.0f} s"))
        T["lens"] = time.time() - t0
        if cand_cache:
            save_candidates(cand_cache, cands, fp)
    n_pairs = cands.groupby(["seg1", "seg2"]).ngroups if len(cands) else 0
    stats.update(candidates=len(cands), pairs=n_pairs)
    log(f"lens: {len(cands)} candidates on {n_pairs} flight pairs" + (f" ({T['lens']:.1f} s)" if "lens" in T else ""))
    if not len(cands):
        out = pd.DataFrame(columns=RESULT_COLUMNS); out.attrs["stats"] = stats
        return out
    w1 = EventWindow.from_fit(sp, ev1, harmonics[0], mode=score_mode)
    w2 = EventWindow.from_fit(sp, ev2, harmonics[1], mode=score_mode) if not single else None
    per_pair = cands.groupby(["seg1", "seg2"]).size().max()
    if pts is None and (n_pairs > exact_pairs or per_pair > exact_points):
        t0 = time.time()
        t1 = build_table(w1, d1, v1, SHIFTS, c, workers, _pool_map)
        t2 = build_table(w2, d2, v2, DELTAS, c, workers, _pool_map) if not single else None
        T["tables"] = time.time() - t0
        t0 = time.time()
        ts = table_scores(cands, t1, t2, params.dt_heard)
        sel = cands.iloc[shortlist(cands, ts, exact_pairs, exact_points)].reset_index(drop=True)
        T["lookup"] = time.time() - t0
        log(f"tables: {len(cands)} candidates looked up ({T['tables']:.1f} s to build); best {exact_points} points of "
            f"the best {exact_pairs} pairs -> {len(sel)} rows on {sel.groupby(['seg1', 'seg2']).ngroups} pairs")
    else:
        sel = cands
    t0 = time.time()
    sc = score_candidates(w1, w2, st, sel, c, workers=workers)
    res = _results(st, sel, sc, "exact")
    T["exact"] = time.time() - t0
    log(f"exact: {len(sel)} rows scored on the stored tracks ({T['exact']:.1f} s)")
    stats["rows_exact"] = len(sel)
    if pts is None and refine_pairs > 0:
        t0 = time.time()
        fine = refine_candidates(st, res, params, elevation, refine_pairs, refine_half_m, refine_step_m)
        if len(fine):
            rf = _results(st, fine, score_candidates(w1, w2, st, fine, c, workers=workers), "refine")
            res = pd.concat([res, rf], ignore_index=True)
            res["gidx"] = np.arange(len(res))
        T["refine"] = time.time() - t0
        stats["rows_refine"] = len(fine)
        log(f"refine: {len(fine)} points on a {refine_step_m:.0f} m lattice around the best {refine_pairs} pairs ({T['refine']:.1f} s)")
    log("stage times: " + ", ".join(f"{k} {v:.1f} s" for k, v in T.items()) + f"; total {sum(T.values()):.1f} s")
    res = res.sort_values("score", ascending=False).reset_index(drop=True)
    res.attrs["stats"] = stats
    return res


# ----------------------------------------------------------------------------------------
# more events: rescore the top pairs with further flyovers locked to the same clock
# ----------------------------------------------------------------------------------------

def pick_flight(st: TrackStore, w: EventWindow, obs_enu, t0_map: float, shift: float, c: float,
                max_range_m: float = 12000.0, exclude=()) -> tuple[int, float]:
    """Which flight does an observer at obs_enu hear as the event in window w, given the
    recording's clock (audio = unix + t0_map) and the time shift chosen for event 1? The
    flights within max_range at that time come from the light index; the best locked score
    wins. Returns (segment index, score), (-1, 0) if nothing is near."""
    t_unix = w.fit.t_cpa - t0_map
    best, best_s = -1, 0.0
    for si in sorted(st.near(obs_enu, t_unix, max_range_m + 150.0 * st.bucket_s) - set(exclude)):
        tr = st.track(si, *w.emission_span(t0_map))
        if tr is None:
            continue
        to, r = modelled_curve(tr, obs_enu, c)
        sc = float(w.sweep_scores(to + t0_map, r, [shift], F0_SCALES).max())
        if sc > best_s:
            best, best_s = si, sc
    return best, best_s


def _rescore_pair(args):
    st: TrackStore = _G["store"]; t1 = _G["ev1_t_cpa"]; xw = _G["xw"]; c = _G["c"]; mr = _G["max_range_m"]
    i0, obs, th1, shift, seg1, seg2 = args
    k = int(np.argmax(_G["score"][i0:i0 + len(obs)]))                      # the pair's best point chooses the flights
    chosen = [pick_flight(st, w, obs[k], t1 - float(th1[k]), float(shift[k]), c, mr, exclude=(seg1, seg2))[0] for w in xw]
    out = np.zeros((len(obs), len(xw)))
    maps = t1 - np.asarray(th1, float)
    for j, (w, si) in enumerate(zip(xw, chosen)):
        if si < 0:
            continue
        spans = [w.emission_span(m) for m in (maps.min(), maps.max())]
        tr = st.track(si, min(s[0] for s in spans), max(s[1] for s in spans))
        if tr is None:
            continue
        for n in range(len(obs)):
            to, r = modelled_curve(tr, obs[n], c)
            out[n, j] = w.sweep_scores(to + maps[n], r, [float(shift[n])], F0_SCALES).max()
    return i0, out, chosen


def rescore_with_extra(sp, ev1, extra, st: TrackStore, res: pd.DataFrame, extra_harmonics=None, top_pairs: int = 200,
                       score_mode="brightness", c: float | None = None, max_range_m: float = 12000.0,
                       workers: int | None = None, log=_log) -> pd.DataFrame:
    """Add further measured events to the top pairs of a search result. For each of the best
    `top_pairs` flight pairs, each extra event is matched to a flight once, at the pair's best
    point (pick_flight), then scored at every point of that pair with that point's own time
    shift, so every event is locked to one clock. Returns those pairs' rows re-ranked by
    score = s1 + s2 + s3 + ..., with the two-event score kept as score12 and the chosen
    flights as x_icao1, x_icao2, ..."""
    c = speed_of_sound() if c is None else c
    workers = default_workers() if workers is None else workers
    extra = list(extra); extra_harmonics = list(extra_harmonics or [1] * len(extra))
    xw = [EventWindow.from_fit(sp, ev, h, mode=score_mode) for ev, h in zip(extra, extra_harmonics)]
    for w in xw:
        w.comb()
    best = res.sort_values("score", ascending=False).drop_duplicates(["seg1", "seg2"]).head(top_pairs)
    keys = set(zip(best.seg1.tolist(), best.seg2.tolist()))
    sel = res[[k in keys for k in zip(res.seg1.tolist(), res.seg2.tolist())]].sort_values(["seg1", "seg2"]).reset_index(drop=True)
    obs = lla_to_enu(sel.lat.to_numpy(), sel.lon.to_numpy(), sel.elev.to_numpy(), st.lat0, st.lon0)
    _G.update(store=st, ev1_t_cpa=float(ev1.t_cpa), xw=xw, c=c, max_range_m=max_range_m, score=sel.score.to_numpy())
    s1a, s2a = sel.seg1.to_numpy(), sel.seg2.to_numpy()
    bounds = np.flatnonzero(np.r_[True, (s1a[1:] != s1a[:-1]) | (s2a[1:] != s2a[:-1]), True])
    jobs = [(int(a), obs[a:b], sel.t_heard1.to_numpy()[a:b], sel["shift"].to_numpy()[a:b], int(s1a[a]), int(s2a[a]))
            for a, b in zip(bounds[:-1], bounds[1:])]
    xs = np.zeros((len(sel), len(xw))); chosen = np.full((len(sel), len(xw)), -1, int)
    t0 = time.time()
    for n, (i0, out, ch) in enumerate(_pool_map(_rescore_pair, jobs, workers, chunksize=1)):
        xs[i0:i0 + len(out)] = out; chosen[i0:i0 + len(out)] = ch
        if (n + 1) % 50 == 0:
            log(f"  rescored {n + 1}/{len(jobs)} pairs ({time.time() - t0:.0f} s)")
    sel["score12"] = sel.score
    for j in range(len(xw)):
        sel[f"s{j + 3}"] = xs[:, j]
        sel[f"x_icao{j + 1}"] = np.where(chosen[:, j] >= 0, st.icao[np.maximum(chosen[:, j], 0)], "")
    sel["score"] = sel.score12 + xs.sum(1)
    log(f"extras: {len(sel)} points on {len(jobs)} pairs rescored with {len(xw)} further events ({time.time() - t0:.1f} s)")
    return sel.sort_values("score", ascending=False).reset_index(drop=True)


# ----------------------------------------------------------------------------------------
# a small lat/lon grid, for refine's fine lattice around one result
# ----------------------------------------------------------------------------------------

@dataclass
class Grid:
    lat: np.ndarray
    lon: np.ndarray
    elev: np.ndarray

    @classmethod
    def make(cls, lat_min, lat_max, lon_min, lon_max, spacing_m, elevation_fn=None):
        lat0 = (lat_min + lat_max) / 2
        dlat = spacing_m / 111320.0
        dlon = spacing_m / (111320.0 * np.cos(np.radians(lat0)))
        la, lo = np.meshgrid(np.arange(lat_min, lat_max, dlat), np.arange(lon_min, lon_max, dlon), indexing="ij")
        la, lo = la.ravel(), lo.ravel()
        el = np.maximum(elevation_fn(la, lo), 0.0) if elevation_fn else np.zeros_like(la)
        return cls(la, lo, el)

    def __len__(self):
        return len(self.lat)

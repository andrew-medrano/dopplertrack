"""The track store: one day of propeller flights in a form a search can read in a second.

Every airborne stretch of every trace (a segment) is kept as knots every KNOT_S seconds, each
with time, position and velocity, and the track between two knots is the cubic Hermite curve
through them. Positions and velocities come from a quadratic least-squares fit over SMOOTH_S
seconds of the 1 Hz track (Savitzky-Golay), so range rates are not differenced from noisy
positions. All positions are in ONE east-north-up frame at (lat0, lon0); lla_to_enu is a true
ECEF rotation, so distances are exact anywhere on Earth.

Next to the knots the store holds the light index the join and the silence filter run on:
every segment sampled once per BUCKET_S seconds at the bucket centre, sorted by bucket.

On disk a store is a directory of .npy files plus meta.json (format name and version, frame,
knot and bucket spacing, sizes). Arrays are memory-mapped on load, so a store opens in well
under a second, costs no RAM until it is read, and is shared by forked worker processes.

    meta.json   format, version, lat0, lon0, t0, b0, knot_s, smooth_s, bucket_s, counts
    kt.npy      (K,)   float32  knot times, s after t0
    kp.npy      (K,3)  float32  knot positions, m
    kv.npy      (K,3)  float32  knot velocities, m/s
    off.npy     (S+1,) int64    segment s owns knots off[s]:off[s+1]
    icao.npy    (S,)   <U6      ICAO address per segment
    trace.npy   (S,)   int32    trace number per segment (several segments can share a flight)
    sb.npy      (M,)   int32    light index: bucket per sample, minus b0, sorted
    ss.npy      (M,)   int32    light index: segment per sample
    spos.npy    (M,3)  float32  light index: position per sample, m
    sspd.npy    (M,)   float32  light index: speed per sample, m/s

A store is built from traces (`TrackStore.build`, what `dopplertrack fetch` does) or from a
1 Hz day index written by dopplertrack 0.3 (`TrackStore.from_day_index`, read through
memory maps so a 6 GB index never has to fit in RAM).
"""
from __future__ import annotations

import json
import os
import struct
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.signal import savgol_filter
from scipy.spatial import cKDTree

from .geo import enu_to_lla

FORMAT = "dopplertrack-store"
VERSION = 1
KNOT_S = 15           # knot spacing, s (samples of the 1 Hz track)
SMOOTH_S = 31         # Savitzky-Golay window, s
BUCKET_S = 10.0       # light index spacing, s
MIN_LEN = 20          # shortest segment kept, s
_ARRAYS = ("kt", "kp", "kv", "off", "icao", "trace", "sb", "ss", "spos", "sspd")


@dataclass
class Track:
    """A stretch of one segment sampled from the store: emission times (unix s), positions and
    velocities (m, m/s, in the store's frame)."""
    t: np.ndarray
    pos: np.ndarray
    vel: np.ndarray

    def __len__(self):
        return len(self.t)

    @property
    def speed(self) -> np.ndarray:
        return np.linalg.norm(self.vel, axis=1)


@dataclass
class Chain:
    """The knots of one segment over a time window, float64: what closest-approach queries run on."""
    t: np.ndarray       # (K,) unix s
    p: np.ndarray       # (K,3)
    v: np.ndarray       # (K,3)

    def __len__(self):
        return len(self.t)


def hermite(t0, t1, p0, p1, v0, v1, t, derivs: int = 1):
    """Cubic Hermite position (and derivatives) at times t inside [t0, t1], vectorised over
    rows: t0, t1, t (n,) and p0, p1, v0, v1 (n, 3). Returns (p, dp/dt) or (p, dp/dt, d2p/dt2)."""
    h = (t1 - t0)[:, None]
    u = ((t - t0)[:, None]) / h
    u2 = u * u; u3 = u2 * u
    p = (2 * u3 - 3 * u2 + 1) * p0 + (u3 - 2 * u2 + u) * h * v0 + (-2 * u3 + 3 * u2) * p1 + (u3 - u2) * h * v1
    dp = ((6 * u2 - 6 * u) * p0 + (3 * u2 - 4 * u + 1) * h * v0 + (-6 * u2 + 6 * u) * p1 + (3 * u2 - 2 * u) * h * v1) / h
    if derivs == 1:
        return p, dp
    ddp = ((12 * u - 6) * p0 + (6 * u - 4) * h * v0 + (-12 * u + 6) * p1 + (6 * u - 2) * h * v1) / (h * h)
    return p, dp, ddp


def segment_knots(t: np.ndarray, enu: np.ndarray, knot_s: int = KNOT_S, smooth_s: int = SMOOTH_S):
    """Knot times, smoothed positions and velocities of one uniformly sampled segment: a
    quadratic least-squares fit over smooth_s samples (or the whole segment when shorter),
    taken every knot_s samples plus the last one."""
    n = len(t)
    dt = float(t[1] - t[0])
    w = min(smooth_s, n if n % 2 else n - 1)
    if w >= 5:
        p = savgol_filter(enu, w, 2, axis=0)
        v = savgol_filter(enu, w, 2, deriv=1, delta=dt, axis=0)
    else:
        p = np.asarray(enu, float); v = np.gradient(enu, t, axis=0)
    k = np.r_[np.arange(0, n - 1, knot_s), n - 1]
    return t[k], p[k], v[k]


class _Builder:
    """Accumulates segments' knots; finish() adds the light index and returns a TrackStore."""

    def __init__(self, lat0: float, lon0: float, knot_s: int = KNOT_S, smooth_s: int = SMOOTH_S):
        self.lat0, self.lon0, self.knot_s, self.smooth_s = float(lat0), float(lon0), knot_s, smooth_s
        self.kt, self.kp, self.kv, self.n, self.icao, self.trace = [], [], [], [], [], []

    def add(self, icao: str, trace: int, t: np.ndarray, enu: np.ndarray, min_len: int = MIN_LEN):
        if len(t) < min_len or not np.all(np.isfinite(enu)):
            return
        kt, kp, kv = segment_knots(t, enu, self.knot_s, self.smooth_s)
        self.kt.append(kt); self.kp.append(kp.astype(np.float32)); self.kv.append(kv.astype(np.float32))
        self.n.append(len(kt)); self.icao.append(icao); self.trace.append(trace)

    def finish(self, bucket_s: float = BUCKET_S) -> "TrackStore":
        kt = np.concatenate(self.kt) if self.kt else np.zeros(0)
        t0 = float(np.floor(kt.min())) if len(kt) else 0.0
        a = dict(kt=(kt - t0).astype(np.float32),
                 kp=np.concatenate(self.kp) if self.kp else np.zeros((0, 3), np.float32),
                 kv=np.concatenate(self.kv) if self.kv else np.zeros((0, 3), np.float32),
                 off=np.r_[0, np.cumsum(self.n)].astype(np.int64),
                 icao=np.array(self.icao, dtype="<U6"), trace=np.array(self.trace, np.int32))
        meta = dict(format=FORMAT, version=VERSION, lat0=self.lat0, lon0=self.lon0, t0=t0, knot_s=self.knot_s,
                    smooth_s=self.smooth_s, bucket_s=float(bucket_s))
        st = TrackStore(meta, a)
        st._add_light_index(bucket_s)
        return st


@dataclass
class TrackStore:
    meta: dict
    a: dict = field(repr=False)
    path: Path | None = None
    _trees: dict = field(default_factory=dict, repr=False)
    _slla: tuple | None = field(default=None, repr=False)

    # ------------------------------------------------------------------ construction
    @classmethod
    def build(cls, traces, lat0: float | None = None, lon0: float | None = None, bucket_s: float = BUCKET_S,
              knot_s: int = KNOT_S, smooth_s: int = SMOOTH_S, min_len: int = MIN_LEN) -> "TrackStore":
        """A store from adsb Traces: every airborne stretch resampled to 1 Hz in the frame at
        (lat0, lon0) (default: the traces' mean position), smoothed, reduced to knots."""
        traces = list(traces)
        if lat0 is None or lon0 is None:
            lat0 = float(np.mean([np.median(tr.lat) for tr in traces])) if traces else 0.0
            lon0 = float(np.mean([np.median(tr.lon) for tr in traces])) if traces else 0.0
        b = _Builder(lat0, lon0, knot_s, smooth_s)
        for i, tr in enumerate(traces):
            for sl in tr.segments():
                t, enu = tr.enu(sl, lat0, lon0, 1.0)
                b.add(tr.icao, i, t, enu, min_len)
        return b.finish(bucket_s)

    @classmethod
    def from_day_index(cls, path, bucket_s: float = BUCKET_S, knot_s: int = KNOT_S, smooth_s: int = SMOOTH_S,
                       progress=None) -> "TrackStore":
        """A store from a 1 Hz day index .npz (dopplertrack 0.3 `fetch` output: lat0, lon0, off, t,
        enu, icao, trace_idx). The big members are memory-mapped, one segment read at a time."""
        z = {k: _npz_member(path, k) for k in ("lat0", "lon0", "off", "t", "enu", "icao", "trace_idx")}
        off = np.asarray(z["off"]); icao = np.asarray(z["icao"]); tix = np.asarray(z["trace_idx"])
        b = _Builder(float(z["lat0"]), float(z["lon0"]), knot_s, smooth_s)
        for s in range(len(off) - 1):
            lo, hi = int(off[s]), int(off[s + 1])
            b.add(str(icao[s]), int(tix[s]), np.asarray(z["t"][lo:hi], float), np.asarray(z["enu"][lo:hi], float), min_len=2)
            if progress and (s + 1) % 5000 == 0:
                progress(s + 1, len(off) - 1)
        return b.finish(bucket_s)

    def _add_light_index(self, bucket_s: float):
        """Every segment sampled at the centre of each bucket_s bucket it spans, sorted by bucket."""
        kt = self.kt_abs(slice(None)); off = self.a["off"]
        S = len(off) - 1
        if S == 0:
            self.meta["b0"] = 0
            self.a.update(sb=np.zeros(0, np.int32), ss=np.zeros(0, np.int32), spos=np.zeros((0, 3), np.float32), sspd=np.zeros(0, np.float32))
            self._count(); return
        first, last = kt[off[:-1]], kt[off[1:] - 1]
        blo = np.ceil((first - bucket_s / 2) / bucket_s).astype(np.int64)
        bhi = np.floor((last - bucket_s / 2) / bucket_s).astype(np.int64)
        n = np.maximum(bhi - blo + 1, 0)
        ss = np.repeat(np.arange(S, dtype=np.int32), n)
        start = np.repeat(np.cumsum(n) - n, n)
        b = np.repeat(blo, n) + (np.arange(len(ss)) - start)
        tq = (b + 0.5) * bucket_s
        # interval of each sample: global knot index i with kt[i] <= tq < kt[i+1] inside its segment
        key = np.arange(S, dtype=np.float64).repeat(np.diff(off)) * 1e7 + (kt - self.meta["t0"])
        i = np.searchsorted(key, ss * 1e7 + (tq - self.meta["t0"]), "right") - 1
        i = np.clip(i, off[ss], off[ss + 1] - 2)
        kp, kv = self.a["kp"], self.a["kv"]
        p, v = hermite(kt[i], kt[i + 1], kp[i].astype(float), kp[i + 1].astype(float), kv[i].astype(float),
                       kv[i + 1].astype(float), tq)
        b0 = int(b.min()) if len(b) else 0
        order = np.argsort(b, kind="stable")
        self.meta["b0"] = b0
        self.a.update(sb=(b[order] - b0).astype(np.int32), ss=ss[order], spos=p[order].astype(np.float32),
                      sspd=np.linalg.norm(v[order], axis=1).astype(np.float32))
        self._count()

    def _count(self):
        self.meta.update(segments=int(len(self.a["off"]) - 1), knots=int(len(self.a["kt"])), samples=int(len(self.a["sb"])))

    # ------------------------------------------------------------------ disk
    def save(self, path) -> Path:
        """Write the store to a directory; meta.json goes last, so a store whose meta.json
        exists is complete."""
        d = Path(path); d.mkdir(parents=True, exist_ok=True)
        (d / "meta.json").unlink(missing_ok=True)
        for k in _ARRAYS:
            np.save(d / f"{k}.npy", np.asarray(self.a[k]))
        tmp = d / f"meta.json.{os.getpid()}.tmp"
        tmp.write_text(json.dumps(self.meta, indent=1))
        os.replace(tmp, d / "meta.json")
        self.path = d
        return d

    @classmethod
    def load(cls, path) -> "TrackStore":
        d = Path(path)
        mp = d / "meta.json"
        if not mp.exists():
            if d.is_file() and d.suffix == ".npz":
                raise ValueError(f"{d} is a 1 Hz day index; convert it with `dopplertrack store {d} STORE_DIR`")
            raise FileNotFoundError(f"{d} is not a track store (no meta.json); write one with `dopplertrack fetch`")
        meta = json.loads(mp.read_text())
        if meta.get("format") != FORMAT:
            raise ValueError(f"{d}: not a dopplertrack track store")
        if meta.get("version") != VERSION:
            raise ValueError(f"{d}: track store format version {meta.get('version')}, this dopplertrack reads "
                             f"version {VERSION}; rebuild it with `dopplertrack fetch` or `dopplertrack store`")
        a = {k: np.load(d / f"{k}.npy", mmap_mode=None if k == "icao" else "r") for k in _ARRAYS}
        st = cls(meta, a, d)
        if (len(a["off"]) - 1, len(a["kt"]), len(a["sb"])) != (meta["segments"], meta["knots"], meta["samples"]):
            raise ValueError(f"{d}: array sizes do not match meta.json; the store is incomplete or damaged, rebuild it")
        return st

    # ------------------------------------------------------------------ basics
    @property
    def lat0(self) -> float:
        return float(self.meta["lat0"])

    @property
    def lon0(self) -> float:
        return float(self.meta["lon0"])

    @property
    def bucket_s(self) -> float:
        return float(self.meta["bucket_s"])

    @property
    def icao(self) -> np.ndarray:
        return self.a["icao"]

    @property
    def sb(self):
        return self.a["sb"]

    @property
    def ss(self):
        return self.a["ss"]

    @property
    def spos(self):
        return self.a["spos"]

    @property
    def sspd(self):
        return self.a["sspd"]

    def __len__(self):
        return int(len(self.a["off"]) - 1)

    def nbytes(self) -> int:
        return int(sum(np.asarray(v).nbytes for v in self.a.values()))

    def fingerprint(self) -> dict:
        m = self.meta
        return {k: m[k] for k in ("version", "lat0", "lon0", "t0", "b0", "knot_s", "bucket_s", "segments", "knots", "samples")}

    def kt_abs(self, sl) -> np.ndarray:
        return np.asarray(self.a["kt"][sl], float) + self.meta["t0"]

    def span(self, s: int) -> tuple[float, float]:
        a, b = int(self.a["off"][s]), int(self.a["off"][s + 1])
        return float(self.a["kt"][a]) + self.meta["t0"], float(self.a["kt"][b - 1]) + self.meta["t0"]

    # ------------------------------------------------------------------ tracks
    def chain(self, s: int, t_a: float = -np.inf, t_b: float = np.inf) -> Chain | None:
        """Knots of segment s covering [t_a, t_b] with one knot of margin each side; None when
        that is fewer than two knots."""
        a, b = int(self.a["off"][s]), int(self.a["off"][s + 1])
        kt = self.a["kt"][a:b]
        t0 = self.meta["t0"]
        i = max(int(np.searchsorted(kt, t_a - t0, "left")) - 1, 0)
        j = min(int(np.searchsorted(kt, t_b - t0, "right")) + 1, b - a)
        if j - i < 2:
            return None
        return Chain(np.asarray(kt[i:j], float) + t0, np.asarray(self.a["kp"][a + i:a + j], float),
                     np.asarray(self.a["kv"][a + i:a + j], float))

    def track(self, s: int, t_a: float = -np.inf, t_b: float = np.inf, dt: float = 1.0) -> Track | None:
        """Segment s sampled every dt seconds (on whole seconds) over [t_a, t_b] ∩ its span;
        None when that is fewer than two samples."""
        ch = self.chain(s, t_a, t_b)
        if ch is None:
            return None
        lo, hi = max(t_a, ch.t[0]), min(t_b, ch.t[-1])
        t = np.arange(np.ceil(lo / dt) * dt, np.floor(hi / dt) * dt + dt / 2, dt)
        if len(t) < 2:
            return None
        p, v = chain_eval(ch, t)
        return Track(t, p, v)

    # ------------------------------------------------------------------ the light index
    def bucket(self, t_unix) -> int:
        """Bucket id (relative to the store's first bucket) holding unix time t."""
        return int(np.floor(t_unix / self.bucket_s)) - int(self.meta["b0"])

    def bucket_start(self, b) -> float:
        return (b + self.meta["b0"]) * self.bucket_s

    def bucket_time(self, b) -> float:
        return (b + self.meta["b0"] + 0.5) * self.bucket_s

    def bucket_slice(self, b: int) -> slice:
        sb = self.a["sb"]
        return slice(int(np.searchsorted(sb, b, "left")), int(np.searchsorted(sb, b, "right")))

    def tree(self, b: int) -> cKDTree | None:
        """KD-tree of every light-index sample in bucket b (cached); None if the bucket is empty."""
        if b not in self._trees:
            sl = self.bucket_slice(b)
            self._trees[b] = cKDTree(np.asarray(self.a["spos"][sl], float)) if sl.stop > sl.start else None
        return self._trees[b]

    def near(self, obs_enu, t_unix: float, radius_m: float) -> set[int]:
        """Segments with a light-index sample within radius of obs_enu in the buckets around t_unix."""
        b = self.bucket(t_unix)
        out: set[int] = set()
        for bb in (b - 1, b, b + 1):
            tr = self.tree(bb)
            if tr is None:
                continue
            sl = self.bucket_slice(bb)
            hit = tr.query_ball_point(np.asarray(obs_enu, float), radius_m)
            out.update(int(x) for x in self.a["ss"][np.asarray(hit, int) + sl.start])
        return out

    def sample_mask_in_bbox(self, bbox, pad_m: float = 0.0) -> np.ndarray:
        """Light-index samples whose position lies inside (lat_min, lat_max, lon_min, lon_max), padded."""
        if self._slla is None:
            p = np.asarray(self.a["spos"], float)
            self._slla = enu_to_lla(p[:, 0], p[:, 1], p[:, 2], self.lat0, self.lon0)[:2]
        la, lo = self._slla
        dla = pad_m / 111320.0; dlo = pad_m / (111320.0 * np.cos(np.radians((bbox[0] + bbox[1]) / 2)))
        return (la >= bbox[0] - dla) & (la <= bbox[1] + dla) & (lo >= bbox[2] - dlo) & (lo <= bbox[3] + dlo)


def chain_eval(ch: Chain, t, derivs: int = 1):
    """Hermite position and velocity (and acceleration with derivs=2) of a chain at times t
    (clamped to its span)."""
    t = np.clip(np.asarray(t, float), ch.t[0], ch.t[-1])
    i = np.clip(np.searchsorted(ch.t, t, "right") - 1, 0, len(ch.t) - 2)
    return hermite(ch.t[i], ch.t[i + 1], ch.p[i], ch.p[i + 1], ch.v[i], ch.v[i + 1], t, derivs)


def _npz_member(path, name: str):
    """One array of an .npz: memory-mapped when the member is stored uncompressed (np.savez),
    read whole otherwise (np.savez_compressed)."""
    with zipfile.ZipFile(path) as zf:
        info = zf.getinfo(name + ".npy")
        stored = info.compress_type == zipfile.ZIP_STORED
    if not stored:
        with np.load(path, allow_pickle=False) as z:
            return z[name]
    with open(path, "rb") as fh:
        fh.seek(info.header_offset)
        h = fh.read(30)
        n_name, n_extra = struct.unpack("<HH", h[26:30])
        fh.seek(info.header_offset + 30 + n_name + n_extra)
        major, _ = np.lib.format.read_magic(fh)
        read = np.lib.format.read_array_header_1_0 if major == 1 else np.lib.format.read_array_header_2_0
        shape, fortran, dtype = read(fh)
        offset = fh.tell()
    if dtype.hasobject:
        raise ValueError(f"{path}:{name} holds Python objects")
    if int(np.prod(shape)) == 0 or shape == ():
        with np.load(path, allow_pickle=False) as z:
            return z[name]
    return np.memmap(path, dtype=dtype, mode="r", shape=shape, order="F" if fortran else "C", offset=offset)

"""High-level steps the CLI and notebooks share: spectrogram, track store, refine."""
from __future__ import annotations

import time

import numpy as np
import pandas as pd

from .adsb import AircraftTypes, load_region, save_traces
from .audio import compute_spectrogram, load_audio
from .cpa import closest_approach
from .doppler import speed_of_sound
from .geo import lla_to_enu
from .score import F0_SCALES, EventWindow, locked_score, modelled_curve
from .search import Grid, best_per_pair, pick_flight
from .store import TrackStore
from .tracks import EventFit

WORLD = (-90.0, 90.0, -180.0, 180.0)

__all__ = ["make_spectrogram", "fetch", "convert_day_index", "best_per_pair", "refine", "log"]


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def make_spectrogram(audio_path, out_path, hop_s=0.25, fmin=40, fmax=450, nperseg=8192, whiten_window_s=600):
    x, sr = load_audio(audio_path)
    sp = compute_spectrogram(x, sr, nperseg=nperseg, hop_s=hop_s, fmin=fmin, fmax=fmax, whiten_window_s=whiten_window_s)
    sp.save(out_path)
    return sp


def _report(st: TrackStore, out_dir, t0):
    log(f"wrote {out_dir}: {len(st)} segments, {st.meta['knots']} knots, {st.meta['samples']} light-index samples, "
        f"{st.nbytes() / 1e6:.0f} MB, origin {st.lat0:.4f}, {st.lon0:.4f} ({time.time() - t0:.0f} s)")


def fetch(date: str, out_dir, bbox=None, archive=None, prop_only=True, pod="planes-readsb-prod-0",
          origin=None, traces_path=None, workers=None) -> TrackStore:
    """One UTC day of adsb.lol history -> the track store a search runs on.

    The traces are streamed from the GitHub release (or read from an extracted archive
    directory), kept if they touch bbox (default: everywhere), stripped of known jets, and
    written as a track store in the frame at `origin` (default: the traces' mean position).
    traces_path also saves the raw traces."""
    bbox = WORLD if bbox is None else tuple(bbox)
    t0 = time.time()
    if archive:
        log(f"scanning {archive} for traces in {bbox}")
        traces = load_region(archive, bbox, workers=workers)
    else:
        from .fetch import fetch_day
        log(f"streaming adsb.lol {date} ({pod}) for traces in {bbox}")
        def prog(n, nbytes):
            log(f"  {n} trace files scanned, {nbytes / 1e9:.2f} GB read, {nbytes / 1e6 / max(time.time() - t0, 1):.1f} MB/s")
        traces = fetch_day(date, bbox, pod=pod, progress=prog)
    log(f"{len(traces)} traces touch the region ({time.time() - t0:.0f} s)")
    if prop_only:
        types = AircraftTypes()
        traces = [t for t in traces if types.is_prop(t.typ)]
        log(f"{len(traces)} after dropping known jets")
    if traces_path:
        save_traces(traces, traces_path)
    t0 = time.time()
    o = origin or (None, None)
    st = TrackStore.build(traces, o[0], o[1])
    st.save(out_dir)
    _report(st, out_dir, t0)
    return st


def convert_day_index(day_npz, out_dir) -> TrackStore:
    """A 1 Hz day index written by dopplertrack 0.3 (`fetch ... day.npz`) -> a track store."""
    t0 = time.time()
    log(f"converting {day_npz}")
    st = TrackStore.from_day_index(day_npz, progress=lambda k, n: log(f"  {k}/{n} segments"))
    st.save(out_dir)
    _report(st, out_dir, t0)
    return st


def refine(sp, ev1: EventFit, ev2: EventFit, st: TrackStore, res_row, spacing_m=50.0, half_m=600.0, harmonics=(1, 1),
           elevation=None, c=None, extra: list[EventFit] | None = None, extra_harmonics: list[int] | None = None,
           score_mode="brightness", max_range_m=12000.0, log=log) -> pd.DataFrame:
    """Fine grid around one result row (a search result's row: seg1, seg2, lat, lon, t_heard1,
    shift) for its flight pair, scored without gates.

    extra: further measured events beyond the two the search used. For each, the flight is
    chosen once at the row's point (pick_flight), then its modelled curve is scored at every
    fine point with the same time shift as event 1, so every event is locked to the same
    clock. A low, lateral pass pins the location far better than a high overhead one, which
    is what the extra events are for. Columns s3, s4, ... are added and summed into score;
    the chosen flights are in the DataFrame's .attrs["extra_icao"]."""
    c = speed_of_sound() if c is None else c
    extra = list(extra or []); extra_harmonics = list(extra_harmonics or [1] * len(extra))
    lat, lon = float(res_row.lat), float(res_row.lon)
    dlat = half_m / 111320.0; dlon = half_m / (111320.0 * np.cos(np.radians(lat)))
    g = Grid.make(lat - dlat, lat + dlat, lon - dlon, lon + dlon, spacing_m, elevation)
    s1, s2 = int(res_row.seg1), int(res_row.seg2)
    w1 = EventWindow.from_fit(sp, ev1, harmonics[0], mode=score_mode)
    w2 = EventWindow.from_fit(sp, ev2, harmonics[1], mode=score_mode)
    el0 = float(elevation(np.array([lat]), np.array([lon]))[0]) if elevation else 0.0
    obs0 = lla_to_enu(lat, lon, el0, st.lat0, st.lon0)
    t0_map0 = ev1.t_cpa - float(res_row.t_heard1)
    xsegs, xwins = [], []
    for ev, h in zip(extra, extra_harmonics):
        w = EventWindow.from_fit(sp, ev, h, mode=score_mode)
        k, sc = pick_flight(st, w, obs0, t0_map0, float(res_row["shift"]), c, max_range_m)
        if k < 0:
            raise ValueError(f"no flight near the candidate at the time of the extra event at {ev.t_cpa:.0f} s")
        log(f"extra event at {ev.t_cpa:.0f} s: {st.icao[k]} (score {sc:.2f} at the candidate)")
        xsegs.append(k); xwins.append(w)
    obs_all = lla_to_enu(g.lat, g.lon, g.elev, st.lat0, st.lon0)
    te1 = float(res_row.t_heard1) - float(res_row.d1) / c
    ch1 = st.chain(s1, te1 - 300.0, te1 + 300.0)
    _, t_heard, _, _ = closest_approach(obs_all, ch1, c)
    rows = []
    for i in range(len(g)):
        obs = obs_all[i]
        t0_map = ev1.t_cpa - float(t_heard[i])
        tr1 = st.track(s1, *w1.emission_span(t0_map)); tr2 = st.track(s2, *w2.emission_span(t0_map))
        tot, sc1, sc2, sh = locked_score(w1, w2, tr1, tr2, obs, float(t_heard[i]), c=c)
        xs = []
        for k, w in zip(xsegs, xwins):
            tr = st.track(k, *w.emission_span(t0_map))
            if tr is None:
                xs.append(0.0); continue
            to, r = modelled_curve(tr, obs, c)
            xs.append(float(w.sweep_scores(to + t0_map, r, [sh], F0_SCALES).max()))
        rows.append((g.lat[i], g.lon[i], tot + sum(xs), sc1, sc2, *xs, sh))
    cols = ["lat", "lon", "score", "s1", "s2"] + [f"s{k + 3}" for k in range(len(extra))] + ["shift"]
    out = pd.DataFrame(rows, columns=cols).sort_values("score", ascending=False).reset_index(drop=True)
    out.attrs["extra_icao"] = [str(st.icao[k]) for k in xsegs]
    return out

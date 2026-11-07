"""The track store: knots every 15 s with smoothed positions and velocities, a 10 s light
index, a versioned on-disk format read through memory maps, and conversion from a 1 Hz day
index of dopplertrack 0.3 without loading it into memory."""
import json

import numpy as np
import pytest

from dopplertrack.adsb import Trace
from dopplertrack.geo import enu_to_lla
from dopplertrack.store import KNOT_S, VERSION, Chain, TrackStore, chain_eval, segment_knots

LAT0, LON0 = 40.8765, -73.525
T0 = 1_700_000_000.0


def _trace(icao, t, enu):
    lat, lon, alt = enu_to_lla(enu[:, 0], enu[:, 1], enu[:, 2], LAT0, LON0)
    return Trace(icao, None, "C172", "synthetic", t, lat, lon, alt, np.full_like(t, 100.0), np.zeros(len(t), bool))


def _turn(v, rate_deg_s, n=900, t0=T0):
    """Constant-speed level flight turning at rate_deg_s, sampled every second."""
    t = np.arange(n, dtype=float)
    w = np.radians(rate_deg_s)
    if w == 0:
        x, y = v * t, 0 * t
    else:
        R = v / w
        x, y = R * np.sin(w * t), R * (1 - np.cos(w * t))
    return t0 + t, np.c_[x - x.mean(), y - y.mean() + 3000.0, np.full(n, 1200.0)]


def _traces():
    return [_trace("aaa001", *_turn(110, 0.0)), _trace("aaa002", *_turn(95, 0.3, t0=T0 + 600)),
            _trace("aaa003", *_turn(60, 3.0, t0=T0 + 200))]


def test_knots_follow_the_track():
    # straight and gently turning flight: the Hermite curve through the knots is the track; a
    # standard-rate turn (3 deg/s) at 60 m/s is smoothed by the 31 s quadratic fit, a few metres
    for v, rate, tol_m, tol_v in ((110, 0.0, 0.01, 0.01), (95, 0.3, 0.5, 0.1), (60, 3.0, 10.0, 5.0)):
        t, enu = _turn(v, rate)
        kt, kp, kv = segment_knots(t, enu)
        assert np.allclose(np.diff(kt)[:-1], KNOT_S)
        p, vel = chain_eval(Chain(kt, kp, kv), t)
        e = np.linalg.norm(p - enu, axis=1)[20:-20]
        assert e.max() < tol_m, (v, rate, e.max())
        assert np.abs(np.linalg.norm(vel, axis=1)[20:-20] - v).max() < tol_v


def test_build_save_load_round_trip(tmp_path):
    st = TrackStore.build(_traces(), LAT0, LON0)
    assert len(st) == 3 and list(st.icao) == ["aaa001", "aaa002", "aaa003"]
    assert st.meta["knots"] == len(st.a["kt"]) < 900 * 3 / 10
    st.save(tmp_path / "s")
    back = TrackStore.load(tmp_path / "s")
    assert back.fingerprint() == st.fingerprint()
    assert isinstance(back.a["kp"], np.memmap)                      # memory-mapped, not read
    for k in ("kt", "kp", "kv", "off", "sb", "ss", "spos", "sspd"):
        assert np.array_equal(np.asarray(back.a[k]), np.asarray(st.a[k])), k
    # a 1 Hz track sampled back from the store matches the input to well under a metre
    t, enu = _turn(110, 0.0)
    tr = back.track(0, T0 + 100, T0 + 200)
    assert len(tr) == 101 and np.allclose(tr.t, T0 + np.arange(100, 201))
    assert np.abs(tr.pos - enu[100:201]).max() < 0.5 and np.allclose(tr.speed, 110, atol=0.05)


def test_light_index_samples_every_segment_at_bucket_centres():
    st = TrackStore.build(_traces(), LAT0, LON0)
    assert np.all(np.diff(st.sb) >= 0)
    for s in range(len(st)):
        rows = np.flatnonzero(st.ss == s)
        tq = np.array([st.bucket_time(b) for b in st.sb[rows]])
        lo, hi = st.span(s)
        assert tq.min() >= lo and tq.max() <= hi and len(rows) >= (hi - lo) // st.bucket_s - 1
        p, v = chain_eval(st.chain(s), tq)
        assert np.abs(st.spos[rows] - p).max() < 0.5
        assert np.abs(st.sspd[rows] - np.linalg.norm(v, axis=1)).max() < 0.01
    s0 = int(st.ss[0])
    assert s0 in st.near(st.spos[0], st.bucket_time(st.sb[0]), 10.0)


def test_chain_and_track_windows():
    st = TrackStore.build(_traces(), LAT0, LON0)
    ch = st.chain(0, T0 + 300, T0 + 400)
    assert ch.t[0] <= T0 + 300 and ch.t[-1] >= T0 + 400 and ch.t[1] > T0 + 300 - 1e-6
    assert st.chain(0, T0 + 5000, T0 + 6000) is None
    assert st.track(0, T0 - 50, T0 + 0.5) is None                  # one sample is not a track
    lo, hi = st.span(0)
    assert st.track(0).t[0] == np.ceil(lo) and st.track(0).t[-1] == np.floor(hi)


def test_load_refuses_other_formats(tmp_path):
    st = TrackStore.build(_traces(), LAT0, LON0)
    d = st.save(tmp_path / "s")
    meta = json.loads((d / "meta.json").read_text())
    (d / "meta.json").write_text(json.dumps({**meta, "version": VERSION + 1}))
    with pytest.raises(ValueError, match="version"):
        TrackStore.load(d)
    (d / "meta.json").write_text(json.dumps({**meta, "knots": meta["knots"] + 1}))
    with pytest.raises(ValueError, match="damaged"):
        TrackStore.load(d)
    with pytest.raises(FileNotFoundError):
        TrackStore.load(tmp_path / "nothing")
    np.savez(tmp_path / "day.npz", x=np.zeros(1))
    with pytest.raises(ValueError, match="dopplertrack store"):
        TrackStore.load(tmp_path / "day.npz")


def _day_index_npz(path, traces, compressed=False):
    """A day index in the dopplertrack 0.3 layout (1 Hz ENU segments plus 10 s samples)."""
    t, enu, off, icao, tix = [], [], [0], [], []
    for i, tr in enumerate(traces):
        for sl in tr.segments():
            tt, e = tr.enu(sl, LAT0, LON0, 1.0)
            t.append(tt); enu.append(e); off.append(off[-1] + len(tt)); icao.append(tr.icao); tix.append(i)
    t = np.concatenate(t); enu = np.concatenate(enu)
    spd = np.linalg.norm(np.gradient(enu, t, axis=0), axis=1)
    save = np.savez_compressed if compressed else np.savez
    save(path, lat0=LAT0, lon0=LON0, bucket_s=10.0, off=np.array(off, np.int64), t=t, enu=enu, spd=spd,
         icao=np.array(icao), trace_idx=np.array(tix, np.int32), sb=np.zeros(0, np.int64), ss=np.zeros(0, np.int32),
         sj=np.zeros(0, np.int32))


@pytest.mark.parametrize("compressed", [False, True])
def test_from_day_index_matches_build(tmp_path, compressed):
    traces = _traces()
    _day_index_npz(tmp_path / "day.npz", traces, compressed)
    conv = TrackStore.from_day_index(tmp_path / "day.npz")
    ref = TrackStore.build(traces, LAT0, LON0)
    assert list(conv.icao) == list(ref.icao) and conv.fingerprint() == ref.fingerprint()
    assert np.allclose(conv.a["kp"], ref.a["kp"], atol=1e-3) and np.allclose(conv.a["kv"], ref.a["kv"], atol=1e-4)

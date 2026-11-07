"""The log-frequency-axis scorer against the linear-axis scorer of dopplertrack 0.3: on the
synthetic recording of test_search, for observers around the truth and every placement of the
sweep, the two must agree (correlation > 0.99) and pick the same best observer and placement."""
import numpy as np
import pytest

from _reference import linear_sweep_scores
from dopplertrack import tracks
from dopplertrack.geo import lla_to_enu
from dopplertrack.score import F0_SCALES, SHIFTS, EventWindow, modelled_curve
from dopplertrack.cpa import closest_approach
from test_search import C, LAT0, LON0, _audio, _index, _seg


@pytest.fixture(scope="module")
def scene():
    st = _index()
    sp = _audio(st, (("aaa001", 120.0), ("aaa002", 130.0)))
    ev1 = tracks.fit_event(sp, 230, 400, 90, 170, c=C, objective="brightness")
    rng = np.random.default_rng(3)
    off = np.r_[[[0.0, 0.0]], rng.uniform(-3000, 3000, (39, 2))]                       # the truth and 39 others
    obs = lla_to_enu(LAT0 + off[:, 1] / 111320.0, LON0 + off[:, 0] / (111320.0 * np.cos(np.radians(LAT0))),
                     np.zeros(len(off)), LAT0, LON0)
    return st, sp, ev1, obs


def _sweeps(st, w, ev1, obs, harmonics_axis):
    s = _seg(st, "aaa001")
    _, th, _, _ = closest_approach(obs, st.chain(s), C)
    fast, ref = [], []
    for o, t in zip(obs, th):
        t0_map = ev1.t_cpa - t
        tr = st.track(s, *w.emission_span(t0_map))
        to, r = modelled_curve(tr, o, C)
        fast.append(w.sweep_scores(to + t0_map, r, SHIFTS, F0_SCALES))
        ref.append(linear_sweep_scores(w, to + t0_map, r, SHIFTS, F0_SCALES))
    return np.array(fast), np.array(ref)


@pytest.mark.parametrize("harmonics", [1, 2])
def test_brightness_agrees_with_the_linear_axis(scene, harmonics):
    st, sp, ev1, obs = scene
    w = EventWindow.from_fit(sp, ev1, harmonics)
    fast, ref = _sweeps(st, w, ev1, obs, harmonics)
    assert np.corrcoef(fast.ravel(), ref.ravel())[0, 1] > 0.99
    assert np.abs(fast - ref).max() < 0.05 * ref.max()
    assert np.argmax(fast.max(axis=(1, 2))) == np.argmax(ref.max(axis=(1, 2))) == 0     # the truth wins under both
    i = np.unravel_index(np.argmax(fast[0]), fast[0].shape); j = np.unravel_index(np.argmax(ref[0]), ref[0].shape)
    assert abs(SHIFTS[i[0]] - SHIFTS[j[0]]) <= 0.5 and abs(F0_SCALES[i[1]] - F0_SCALES[j[1]]) <= 0.005


def test_ridge_mode_agrees_with_the_linear_axis(scene):
    st, sp, ev1, obs = scene
    w = EventWindow.from_fit(sp, ev1, 1, mode="ridge")
    fast, ref = _sweeps(st, w, ev1, obs, 1)
    assert np.corrcoef(fast.ravel(), ref.ravel())[0, 1] > 0.99
    assert np.argmax(fast.max(axis=(1, 2))) == np.argmax(ref.max(axis=(1, 2))) == 0

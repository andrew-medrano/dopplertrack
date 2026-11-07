"""Closest approach to a stored track must be the closest approach to the stored curve, on
straight passes and on tracks that turn, circle and come back: checked against brute force on
the Hermite curve sampled every 0.05 s."""
import numpy as np

from dopplertrack.cpa import closest_approach, closest_approach_dense
from dopplertrack.store import Chain, chain_eval, segment_knots

C = 343.0


def _chain(heading_rate, v=100.0, n=600, t0=1.7e9):
    """A level track at speed v whose heading changes at heading_rate (rad/s, one per second)."""
    t = t0 + np.arange(n, dtype=float)
    hd = np.cumsum(heading_rate)
    x = np.cumsum(np.cos(hd) * v); y = np.cumsum(np.sin(hd) * v)
    return Chain(*segment_knots(t, np.c_[x, y, np.full(n, 1000.0)])), x, y


def _points(x, y, n, rng, pad=12000.0):
    return np.c_[rng.uniform(x.min() - pad, x.max() + pad, n), rng.uniform(y.min() - pad, y.max() + pad, n), np.zeros(n)]


def _check(G, ch):
    d, th, sp, inside = closest_approach(G, ch, C)
    D, TH, SP, INSIDE = closest_approach_dense(G, ch, C)
    # never farther than the brute force (which samples the curve every 0.05 s, ~5 m)
    assert (d - D).max() < 0.3, (d - D).max()
    assert (D - d).max() < 1.0
    # and the reported time is where the curve really is that close
    p, _ = chain_eval(ch, th - d / C)
    assert np.abs(np.linalg.norm(p - G, axis=1) - d).max() < 1e-6
    return d, th, sp, inside, D, TH, SP, INSIDE


def test_straight_pass_matches_brute_force():
    rng = np.random.default_rng(0)
    ch, x, y = _chain(np.zeros(600))
    G = _points(x, y, 3000, rng)
    d, th, sp, inside, D, TH, SP, INSIDE = _check(G, ch)
    m = INSIDE & inside
    assert np.abs(th - TH)[m].max() < 0.06 and np.abs(sp - 100.0).max() < 1e-6
    assert (inside == INSIDE).mean() > 0.999


def test_turning_and_circling_tracks_match_brute_force():
    # random heading changes, and every third track circling: other parts of the track come back
    # near most points, which a nearest-knot search alone gets wrong by hundreds of metres
    rng = np.random.default_rng(1)
    for trial in range(30):
        rate = np.cumsum(rng.normal(0, 0.003, 600)) + (0.02 if trial % 3 == 0 else 0.0)
        ch, x, y = _chain(rate, v=60 + 60 * rng.random())
        _check(_points(x, y, 2000, rng), ch)


def test_chord_mode_is_within_the_arc_bound():
    rng = np.random.default_rng(2)
    ch, x, y = _chain(np.full(600, 0.02))
    G = _points(x, y, 2000, rng)
    d, *_ = closest_approach(G, ch, C)
    dc, *_ = closest_approach(G, ch, C, exact=False)
    h = np.diff(ch.t)[:, None]
    D = ch.p[1:] - ch.p[:-1]
    E = 4 / 27 * (np.linalg.norm(h * ch.v[:-1] - D, axis=1) + np.linalg.norm(h * ch.v[1:] - D, axis=1))
    assert np.abs(dc - d).max() <= E.max() + 1e-6


def test_points_past_the_window_ends_are_not_interior():
    ch, x, y = _chain(np.zeros(600))
    G = np.array([[x[0] - 5000, 0, 0], [x[300], 2000, 0], [x[-1] + 5000, 0, 0]])
    d, th, sp, inside = closest_approach(G, ch, C)
    assert list(inside) == [False, True, False]
    assert abs(d[1] - np.hypot(2000, 1000)) < 0.5

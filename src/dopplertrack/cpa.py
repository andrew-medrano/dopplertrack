"""Closest point of approach of many ground points to one stored track.

A track window is a chain of Hermite knots (store.Chain). The closest approach is found in two
steps, both vectorised over the points:

1. Chords. The straight chords between consecutive knots are a polyline; for each point the
   exact nearest chord is found among the CHORD_K chords with the nearest midpoints. That
   candidate list is certified: no chord whose midpoint is farther than the K-th nearest one
   can be closer than (that distance - the longest half-chord). Points the certificate does
   not cover (where another part of a curving or circling track comes back near them) are
   solved against every chord of the window.
2. Curve. A chord and its Hermite arc differ by at most E = 4/27 (|h v0 - D| + |h v1 - D|)
   (h the knot spacing, D the chord vector). Every chord whose distance minus E is below the
   best arc distance found so far is refined on its arc (safeguarded Newton steps on the
   squared distance), so the answer is the closest approach to the stored curve itself, not
   to its polyline, wherever the track turns.

closest_approach_dense is the brute-force reference: the curve sampled every 0.05 s.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from .store import Chain, chain_eval, hermite

CHORD_K = 6          # nearest-midpoint chords examined per point before the certificate
NEWTON_STEPS = 4


def _chords(ch: Chain):
    P0 = ch.p[:-1]; D = ch.p[1:] - ch.p[:-1]
    h = np.diff(ch.t)
    L2 = np.maximum((D * D).sum(1), 1e-9)
    E = 4.0 / 27.0 * (np.linalg.norm(h[:, None] * ch.v[:-1] - D, axis=1) + np.linalg.norm(h[:, None] * ch.v[1:] - D, axis=1))
    return P0, D, L2, E


def _chord_dist(G, P0, D, L2, k):
    """Distance from points G (n,3) to chords k (n, m) and the chord parameter s in [0, 1]."""
    rel = G[:, None, :] - P0[k]
    Dk = D[k]
    s = np.clip((rel * Dk).sum(2) / L2[k], 0.0, 1.0)
    q = rel - s[..., None] * Dk
    return np.sqrt((q * q).sum(2)), s


def _arc_min(G, ch: Chain, k, s):
    """Refine each point's closest approach on the arc of interval k (n,), starting at chord
    parameter s: safeguarded Newton on |p(t) - G|^2 inside [t_k, t_k+1]. Returns (d, t)."""
    t0, t1 = ch.t[k], ch.t[k + 1]
    p0, p1, v0, v1 = ch.p[k], ch.p[k + 1], ch.v[k], ch.v[k + 1]
    t = t0 + s * (t1 - t0)
    p, _ = hermite(t0, t1, p0, p1, v0, v1, t)
    best_t = t; best_d = np.linalg.norm(p - G, axis=1)
    for _ in range(NEWTON_STEPS):
        p, dp, ddp = hermite(t0, t1, p0, p1, v0, v1, t, derivs=2)
        r = p - G
        g = (r * dp).sum(1)
        H = (dp * dp).sum(1) + (r * ddp).sum(1)
        step = np.where(H > 1e-9, -g / np.where(H > 1e-9, H, 1.0), -np.sign(g) * 0.25 * (t1 - t0))
        t = np.clip(t + step, t0, t1)
        p, _ = hermite(t0, t1, p0, p1, v0, v1, t)
        d = np.linalg.norm(p - G, axis=1)
        m = d < best_d
        best_d = np.where(m, d, best_d); best_t = np.where(m, t, best_t)
    return best_d, best_t


def closest_approach(G, ch: Chain, c: float, exact: bool = True):
    """(d, t_heard, speed, interior) of each point of G (n, 3) to the chain's track.

    d is the slant distance at closest approach, t_heard the emission time of closest approach
    plus d / c, speed the track speed there. interior is False where the closest approach is
    at either end of the window, i.e. the true one lies outside it. exact=False stops after
    the chord step: distances are then those to the polyline, within E of the curve's (well
    under a metre on straight flight, tens of metres in a tight turn), which is what a coarse
    gate with slack needs."""
    G = np.asarray(G, float)
    n = len(G)
    if n == 0:
        z = np.zeros(0)
        return z, z, z, np.zeros(0, bool)
    P0, D, L2, E = _chords(ch)
    nk = len(P0)
    kq = min(nk, CHORD_K)
    if kq == nk:
        kk = np.broadcast_to(np.arange(nk)[None, :], (n, nk))
        far = np.full(n, np.inf)                       # every chord examined
    else:
        dm, kk = cKDTree(P0 + 0.5 * D).query(G, k=kq, workers=1)
        far = dm[:, -1] - 0.5 * float(np.sqrt(L2.max()))   # lower bound on the distance to any chord not examined
    dc, sc = _chord_dist(G, P0, D, L2, kk)            # (n, kq)
    j = np.argmin(dc, 1); r = np.arange(n)
    bd, bk, bs = dc[r, j], np.asarray(kk)[r, j].copy(), sc[r, j]
    if not exact:
        bad = np.flatnonzero(far < bd)
        if len(bad):
            d_all, s_all = _chord_dist(G[bad], P0, D, L2, np.broadcast_to(np.arange(nk)[None, :], (len(bad), nk)))
            jj = np.argmin(d_all, 1); rr = np.arange(len(bad))
            bd[bad], bk[bad], bs[bad] = d_all[rr, jj], jj, s_all[rr, jj]
        h = ch.t[bk + 1] - ch.t[bk]
        interior = ~(((bk == 0) & (bs <= 0.0)) | ((bk == nk - 1) & (bs >= 1.0)))
        return bd, ch.t[bk] + bs * h + bd / c, np.sqrt(L2[bk]) / h, interior
    # the arc of the nearest chord is an upper bound; any chord whose distance minus its arc
    # bound E is below it is a contender and is refined on its own arc
    d_best, t_best = _arc_min(G, ch, bk, bs)
    bad = far - float(E.max()) < d_best               # the examined chords may not hold the answer
    blocks = []
    good = np.flatnonzero(~bad)
    if len(good):
        blocks.append((good, dc[good], np.asarray(kk)[good], sc[good]))
    badi = np.flatnonzero(bad)
    if len(badi):
        allk = np.broadcast_to(np.arange(nk)[None, :], (len(badi), nk))
        d_all, s_all = _chord_dist(G[badi], P0, D, L2, allk)
        blocks.append((badi, d_all, np.asarray(allk), s_all))
    for rows, dcm, kkm, scm in blocks:
        cand = (dcm - E[kkm] < d_best[rows][:, None]) & (kkm != bk[rows][:, None])
        if not cand.any():
            continue
        ri, ci = np.nonzero(cand)
        pts = rows[ri]
        d2, t2 = _arc_min(G[pts], ch, kkm[ri, ci], scm[ri, ci])
        order = np.lexsort((d2, pts))                 # several contenders per point: keep the closest
        pts, d2, t2 = pts[order], d2[order], t2[order]
        first = np.r_[True, pts[1:] != pts[:-1]]
        pts, d2, t2 = pts[first], d2[first], t2[first]
        m = d2 < d_best[pts]
        d_best[pts[m]] = d2[m]; t_best[pts[m]] = t2[m]
    _, vel = chain_eval(ch, t_best)
    eps = 1e-6
    interior = (t_best > ch.t[0] + eps) & (t_best < ch.t[-1] - eps)
    return d_best, t_best + d_best / c, np.linalg.norm(vel, axis=1), interior


def closest_approach_dense(G, ch: Chain, c: float, step_s: float = 0.05):
    """Reference for closest_approach: the Hermite curve sampled every step_s seconds and the
    nearest sample taken by brute force."""
    G = np.asarray(G, float)
    t = np.arange(ch.t[0], ch.t[-1] + step_s / 2, step_s)
    t[-1] = min(t[-1], ch.t[-1])
    p, v = chain_eval(ch, t)
    best_d = np.full(len(G), np.inf); best_j = np.zeros(len(G), int)
    for a in range(0, len(t), 2000):
        d = np.linalg.norm(G[:, None, :] - p[None, a:a + 2000], axis=2)
        j = np.argmin(d, 1); dj = d[np.arange(len(G)), j]
        m = dj < best_d
        best_d[m] = dj[m]; best_j[m] = j[m] + a
    tb = t[best_j]
    interior = (best_j > 0) & (best_j < len(t) - 1)
    return best_d, tb + best_d / c, np.linalg.norm(v[best_j], axis=1), interior

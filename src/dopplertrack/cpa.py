"""Closest point of approach of many ground points to one stored track.

A track window is a chain of Hermite knots (store.Chain). The closest approach is found in two
steps, both vectorised over the points:

1. Chords. The straight chords between consecutive knots are a polyline, and the nearest
   chord to each point is found exactly: against every chord (as matrix products) for a
   window of up to DENSE_CHORDS chords, and for longer windows among the chords with the
   nearest midpoints, with a certificate that no other chord can be nearer, falling back to
   every chord for the points it does not cover. A curving or circling track that comes back
   near a point is handled like any other: there is no "nearest knot" shortcut to miss the
   other pass.
2. Curve. A chord and its Hermite arc differ by at most E = 4/27 (|h v0 - D| + |h v1 - D|)
   (h the knot spacing, D the chord vector). The nearest chord's arc, and every other chord
   whose distance minus its E is below the best arc distance found, is refined on the arc
   (safeguarded Newton steps on the squared distance), so the answer is the closest approach
   to the stored curve itself, not to its polyline, wherever the track turns.

closest_approach_dense is the brute-force reference: the curve sampled every 0.05 s.
"""
from __future__ import annotations

import numpy as np

from scipy.spatial import cKDTree

from .store import Chain, chain_eval

NEWTON_STEPS = 4
DENSE_CHORDS = 40         # chains up to this many chords: every point against every chord
NEAR_CHORDS = 8           # longer chains: the chords with the nearest midpoints, then the certificate
_BLOCK = 1 << 17          # points x chords per block of the chord table, to stay in cache


def _chords(ch: Chain):
    P0 = ch.p[:-1]; D = ch.p[1:] - ch.p[:-1]
    h = np.diff(ch.t)
    L2 = np.maximum((D * D).sum(1), 1e-9)
    E = 4.0 / 27.0 * (np.linalg.norm(h[:, None] * ch.v[:-1] - D, axis=1) + np.linalg.norm(h[:, None] * ch.v[1:] - D, axis=1))
    return P0, D, L2, E


def arc_bound(ch: Chain) -> float:
    """Largest distance between the stored curve and its chords over the chain, m."""
    return float(_chords(ch)[3].max())


def _chord_table(G, P0, D, L2):
    """(d, s): distance from every point G (n,3) to every chord (n, nk), and the chord
    parameter of the nearest point in [0, 1]. |G - P0 - sD|^2 expanded into matrix products,
    in coordinates centred on the chain so nothing large cancels."""
    ctr = P0.mean(0)
    Gc = G - ctr; Pc = P0 - ctr
    dot = Gc @ D.T - (Pc * D).sum(1)[None, :]                                  # (G - P0) . D
    s = np.clip(dot / L2[None, :], 0.0, 1.0)
    r2 = (Gc * Gc).sum(1)[:, None] - 2.0 * (Gc @ Pc.T) + (Pc * Pc).sum(1)[None, :]
    d2 = r2 - 2.0 * s * dot + s * s * L2[None, :]
    return np.sqrt(np.maximum(d2, 0.0)), s


def _arc_min(G, ch: Chain, k, s):
    """Refine each point's closest approach on the arc of interval k (n,), starting at chord
    parameter s: safeguarded Newton steps on |p(t) - G|^2 inside [t_k, t_k+1], on the arc's
    power-basis cubic p = a + b u + c u^2 + d u^3 (u = t - t_k). Returns (d, t) at the closest
    place visited."""
    h = (ch.t[k + 1] - ch.t[k])
    hc = h[:, None]
    p0, p1, v0, v1 = ch.p[k], ch.p[k + 1], ch.v[k], ch.v[k + 1]
    A = p0 - G
    B = v0
    Cc = (3.0 * (p1 - p0) / hc - 2.0 * v0 - v1) / hc
    Dd = (2.0 * (p0 - p1) / hc + v0 + v1) / (hc * hc)
    u = s * h
    best_d = np.full(len(G), np.inf); best_u = u
    for it in range(NEWTON_STEPS + 1):
        uc = u[:, None]
        r = A + uc * (B + uc * (Cc + uc * Dd))
        d = np.sqrt((r * r).sum(1))
        m = d < best_d
        best_d = np.where(m, d, best_d); best_u = np.where(m, u, best_u)
        if it == NEWTON_STEPS:
            break
        dp = B + uc * (2.0 * Cc + 3.0 * uc * Dd)
        ddp = 2.0 * Cc + 6.0 * uc * Dd
        g = (r * dp).sum(1)
        H = (dp * dp).sum(1) + (r * ddp).sum(1)
        step = np.where(H > 1e-9, -g / np.where(H > 1e-9, H, 1.0), -np.sign(g) * 0.25 * h)
        u_new = np.clip(u + step, 0.0, h)
        if np.abs(u_new - u).max() < 1e-4:
            u = u_new
            if it < NEWTON_STEPS - 1:      # converged: one last evaluation at the new place
                uc = u[:, None]; r = A + uc * (B + uc * (Cc + uc * Dd)); d = np.sqrt((r * r).sum(1))
                m = d < best_d; best_d = np.where(m, d, best_d); best_u = np.where(m, u, best_u)
            break
        u = u_new
    return best_d, ch.t[k] + best_u


def _chord_gather(G, P0, D, L2, kk):
    """(d, s) from points G (n,3) to the chords kk (n, m) of each point."""
    rel = G[:, None, :] - P0[kk]
    Dk = D[kk]
    s = np.clip((rel * Dk).sum(2) / L2[kk], 0.0, 1.0)
    q = rel - s[..., None] * Dk
    return np.sqrt((q * q).sum(2)), s


def _solve(G, ch: Chain, c: float, exact: bool, kk, dc, sc, L2, E):
    """Closest approach given each point's candidate chords kk (n, m) with their distances dc
    and parameters sc; the answer is among them (or within their arc bounds)."""
    n, nk = len(G), len(ch.t) - 1
    r = np.arange(n)
    j = np.argmin(dc, 1)
    bk = kk[r, j]; bd = dc[r, j]; bs = sc[r, j]
    if not exact:
        t = ch.t[bk] + bs * (ch.t[bk + 1] - ch.t[bk])
        sp = np.linalg.norm(ch.v, axis=1)                    # speed interpolated between knots, so chords
        speed = sp[bk] + bs * (sp[bk + 1] - sp[bk])          # that tie at a shared knot agree
        interior = ~(((bk == 0) & (bs <= 0.0)) | ((bk == nk - 1) & (bs >= 1.0)))
        return bd, t + bd / c, speed, interior
    d_best, t_best = _arc_min(G, ch, bk, bs)
    cand = (dc - E[kk] < d_best[:, None]) & (kk != bk[:, None])
    if cand.any():
        pts, ci = np.nonzero(cand)
        d2, t2 = _arc_min(G[pts], ch, kk[pts, ci], sc[pts, ci])
        order = np.lexsort((d2, pts))                        # several contenders per point: keep the closest
        pts, d2, t2 = pts[order], d2[order], t2[order]
        first = np.r_[True, pts[1:] != pts[:-1]]
        pts, d2, t2 = pts[first], d2[first], t2[first]
        m = d2 < d_best[pts]
        d_best[pts[m]] = d2[m]; t_best[pts[m]] = t2[m]
    _, vel = chain_eval(ch, t_best)
    eps = 1e-6
    interior = (t_best > ch.t[0] + eps) & (t_best < ch.t[-1] - eps)
    return d_best, t_best + d_best / c, np.linalg.norm(vel, axis=1), interior


def _dense(G, ch, c, exact, P0, D, L2, E):
    """Every point against every chord, in blocks that stay in cache."""
    nk = len(P0)
    step = max(_BLOCK // nk, 256)
    out = []
    for a in range(0, len(G), step):
        g = G[a:a + step]
        dc, sc = _chord_table(g, P0, D, L2)
        out.append(_solve(g, ch, c, exact, np.broadcast_to(np.arange(nk)[None, :], dc.shape), dc, sc, L2, E))
    return out[0] if len(out) == 1 else tuple(np.concatenate(x) for x in zip(*out))


def closest_approach(G, ch: Chain, c: float, exact: bool = True):
    """(d, t_heard, speed, interior) of each point of G (n, 3) to the chain's track.

    d is the slant distance at closest approach, t_heard the emission time of closest approach
    plus d / c, speed the track speed there. interior is False where the closest approach is
    at either end of the window, i.e. the true one lies outside it. exact=False stops after
    the chord step: distances are then those to the polyline, within E of the curve's (well
    under a metre on straight flight, tens of metres in a tight turn), which is what a coarse
    gate with slack needs.

    Chains of up to DENSE_CHORDS chords test every chord. Longer ones (two flights together
    for an hour) test the NEAR_CHORDS chords with the nearest midpoints, which provably hold
    the answer when the farthest of them, less the longest half-chord (and twice the largest
    arc bound), is still farther than the nearest; the points where that fails, where the
    track comes back near them, test every chord."""
    G = np.asarray(G, float)
    if len(G) == 0:
        z = np.zeros(0)
        return z, z, z, np.zeros(0, bool)
    P0, D, L2, E = _chords(ch)
    nk = len(P0)
    if nk <= DENSE_CHORDS:
        return _dense(G, ch, c, exact, P0, D, L2, E)
    dm, kk = cKDTree(P0 + 0.5 * D).query(G, k=NEAR_CHORDS, workers=1)
    dc, sc = _chord_gather(G, P0, D, L2, kk)
    far = dm[:, -1] - 0.5 * float(np.sqrt(L2.max())) - (2.0 * float(E.max()) if exact else 0.0)
    bad = far < dc.min(1)
    if not bad.any():
        return _solve(G, ch, c, exact, kk, dc, sc, L2, E)
    out = [np.empty(len(G)), np.empty(len(G)), np.empty(len(G)), np.empty(len(G), bool)]
    good = ~bad
    if good.any():
        for o, x in zip(out, _solve(G[good], ch, c, exact, kk[good], dc[good], sc[good], L2, E)):
            o[good] = x
    for o, x in zip(out, _dense(G[bad], ch, c, exact, P0, D, L2, E)):
        o[bad] = x
    return tuple(out)


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

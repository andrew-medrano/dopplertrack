"""Straight-pass score tables: a cheap first score for every candidate.

Over the few seconds that decide a Doppler curve most flights fly nearly straight and level, so
the curve an observer hears is close to the ideal straight pass with the candidate's CPA
distance and speed. For each event the score of that ideal curve is tabulated once over a
log-spaced grid of distance x speed x time shift (heard CPA pinned at the event's measured
t_cpa, best over the f0 scales), and every candidate gets its score by lookup. The search then
scores exactly, on the real trajectories, only the candidates the tables rank highest.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .doppler import ideal_curve
from .score import F0_SCALES, EventWindow

STEP = 0.015                            # grid step in log(distance) and log(speed): 1.5 %
DELTAS = np.arange(-16, 16.01, 0.5)     # event-2 offsets: geometry gap error plus event-1 shift, s


@dataclass
class ScoreTable:
    d: np.ndarray          # (D,) CPA slant distances, m (log-spaced)
    v: np.ndarray          # (V,) speeds, m/s (log-spaced)
    shifts: np.ndarray     # (S,) time shifts, s
    s: np.ndarray          # (D, V, S) score of the straight pass, best over f0 scales

    def _i(self, g, x):
        return np.clip(np.rint((np.log(np.maximum(x, 1e-9)) - np.log(g[0])) / STEP).astype(np.int64), 0, len(g) - 1)

    def lookup(self, d, v) -> np.ndarray:
        """(n, S) scores at every shift for candidates with CPA distances d and speeds v."""
        return self.s[self._i(self.d, d), self._i(self.v, v)]

    def lookup_at(self, d, v, shift) -> np.ndarray:
        """(n,) score at each candidate's own shift (nearest grid shift, clamped)."""
        sh = np.asarray(shift, float)
        step = self.shifts[1] - self.shifts[0]
        k = np.clip(np.rint((sh - self.shifts[0]) / step).astype(np.int64), 0, len(self.shifts) - 1)
        return self.s[self._i(self.d, d), self._i(self.v, v), k]


_T: dict = {}


def _row(i):
    w, d, vg, sh, c, t = _T["w"], _T["d"][i], _T["v"], _T["shifts"], _T["c"], _T["t"]
    out = np.zeros((len(vg), len(sh)), np.float32)
    for j, v in enumerate(vg):
        r = ideal_curve(t, 1.0, v, d, w.fit.t_cpa - d / c, c=c)     # heard CPA at the event's t_cpa
        out[j] = w.sweep_scores(t, r, sh, F0_SCALES).max(1)
    return i, out


def build_table(w: EventWindow, d_lim, v_lim, shifts, c: float, workers: int = 1, pool_map=None) -> ScoreTable:
    """The score table of one event window over distances d_lim and speeds v_lim (each
    widened by 3 %) and the given shifts. pool_map(fn, items, workers) runs rows in parallel."""
    dg = np.exp(np.arange(np.log(d_lim[0] * 0.97), np.log(d_lim[1] * 1.03), STEP))
    vg = np.exp(np.arange(np.log(max(v_lim[0], 1.0) * 0.97), np.log(v_lim[1] * 1.03), STEP))
    shifts = np.asarray(shifts, float)
    pad = float(np.abs(shifts).max()) + 30.0
    w.comb()                                                          # built once, before any fork
    _T.update(w=w, d=dg, v=vg, shifts=shifts, c=c, t=np.arange(w.fit.t_start - pad, w.fit.t_end + pad, 0.5))
    S = np.zeros((len(dg), len(vg), len(shifts)), np.float32)
    rows = pool_map(_row, range(len(dg)), workers) if pool_map else map(_row, range(len(dg)))
    for i, row in rows:
        S[i] = row
    return ScoreTable(dg, vg, shifts, S)

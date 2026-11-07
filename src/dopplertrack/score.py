"""Score candidate observers by modelling the exact Doppler curve they would hear.

The scorer works on a log-frequency axis (EventWindow.comb): the window's evidence is resampled
once onto (query time x log frequency) with every harmonic summed, so an f0 scale is an index
offset and each (time shift, f0 scale) placement of a modelled curve is one integer gather.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .audio import Spectrogram
from .doppler import speed_of_sound
from .store import Track
from .tracks import EventFit, local_whiten, ridge

F0_SCALES = np.arange(0.95, 1.0501, 0.005)   # f0 scales every placement is swept over
# event-1 time shifts, s. Event 2 is placed with event 1's best shift, so the grid step is also how
# far event 2 can sit from its best placement: at 0.5 s a 0.1 s change in a candidate's heard CPA
# time could move event 1's pick by a step and event 2's score by half a point.
SHIFTS = np.arange(-8, 8.01, 0.25)
DLOG = 0.001                                 # log-frequency axis step (0.1 % of frequency)
_PAD = 64                                    # zero margin of the axis, columns (f0 scales up to 6.6 %)
RIDGE_MIN_DB = 3.0          # ridge frames dimmer than this (whitened dB) carry no evidence
RIDGE_MIN_UNIQUE_DB = 3.0   # ... and so do frames where another path through the window is nearly as bright
RIDGE_EXCLUDE_BINS = 3      # half-width of the band blanked around the ridge to find that other path
RIDGE_SIGMA = 0.015         # frequency tolerance of the ridge score, as a fraction of the ridge frequency
SCORE_MODES = ("brightness", "ridge")


@dataclass
class EventWindow:
    """Whitened spectrogram crop for one event plus its measured parameters.

    mode="brightness" scores a modelled curve by the brightness under it. mode="ridge" first
    extracts the ridge of the measured curve (a Viterbi path over the fit window; harmonics
    are multiples of it) and scores a modelled curve by how closely it follows that ridge:
    each reliable frame contributes exp(-0.5 (df / (RIDGE_SIGMA f))^2), so the score is the
    fraction of the window where the model runs along the ridge, in [0, 1]. A frame is
    reliable when the ridge is at least RIDGE_MIN_DB bright and at least RIDGE_MIN_UNIQUE_DB
    brighter than the best path that avoids it: a bird burst or a second tone offers several
    equally good paths and says nothing about where the aircraft's tone is. That is what a
    faint curve needs: brightness cannot separate candidates on a tone a few dB above the
    noise, but the ridge's position can."""
    fit: EventFit
    sp: Spectrogram
    harmonics: int = 1
    mode: str = "brightness"
    ridge_sigma: float = RIDGE_SIGMA
    ridge_data: tuple | None = field(default=None, repr=False, compare=False)   # (t, f[h], w[h])
    _comb: dict = field(default_factory=dict, repr=False, compare=False)       # n -> comb(n)

    @classmethod
    def from_fit(cls, sp: Spectrogram, fit: EventFit, harmonics: int = 1, pad_s: float = 40.0,
                 mode: str = "brightness", ridge_min_db: float = RIDGE_MIN_DB,
                 ridge_min_unique_db: float = RIDGE_MIN_UNIQUE_DB, ridge_sigma: float = RIDGE_SIGMA):
        if mode not in SCORE_MODES:
            raise ValueError(f"unknown score mode {mode!r}; use one of {SCORE_MODES}")
        w = sp.crop(fit.t_start - pad_s, fit.t_end + pad_s, max(fit.f_post * 0.8, sp.f[0]),
                    min(fit.f_pre * 1.25 * harmonics, sp.f[-1]))
        ew = cls(fit, local_whiten(w), harmonics, mode, ridge_sigma)
        if mode == "ridge":
            ew.extract_ridge(ridge_min_db, ridge_min_unique_db)
        return ew

    def extract_ridge(self, min_db: float = RIDGE_MIN_DB, min_unique_db: float = RIDGE_MIN_UNIQUE_DB):
        """Ridge of the fundamental over the fit window, within the band the fit's asymptotes
        bound (with margin). Harmonic h is h times that ridge (the harmonics of one pass are
        one event; tracking them separately lets a harmonic band latch onto the brighter
        fundamental). A frame is reliable when the ridge is at least min_db bright and at
        least min_unique_db brighter than the best path through the window that stays
        RIDGE_EXCLUDE_BINS bins away from it; harmonic h additionally needs h * ridge to be
        at least min_db bright."""
        w = self.sp.crop(self.fit.t_start, self.fit.t_end, max(self.fit.f_post * 0.8, self.sp.f[0]),
                         min(self.fit.f_pre * 1.25, self.sp.f[-1]))
        fr = ridge(w)
        j = np.rint((fr - w.f[0]) / w.df).astype(int)
        S2 = np.where(np.abs(np.arange(len(w.f))[:, None] - j[None, :]) <= RIDGE_EXCLUDE_BINS, 0.0, w.S)
        fr2 = ridge(Spectrogram(w.t, w.f, S2.astype(np.float32)))
        unique = (w.sample(w.t, fr) - w.sample(w.t, fr2)) >= min_unique_db
        fs = np.array([h * fr for h in range(1, self.harmonics + 1)])
        ws = np.array([(unique & (self.sp.sample(w.t, h * fr) >= min_db)).astype(float)
                       for h in range(1, self.harmonics + 1)])
        self.ridge_data = (w.t, fs, ws)

    def comb(self, n: int = 160):
        """(tq, C, lf0, L): the evidence of the window on n query times of the fit window and a
        log-frequency axis lf0 + DLOG * i (i < L) of the *fundamental*, with every harmonic
        already summed: C[k, i] = sum_h evidence_h(tq[k], h * exp(lf0 + DLOG * i)), stored
        flattened with zero margins so that placements off the axis read zero. Built once per
        window and n.

        On a log axis an f0 scale is a constant index offset, so each (shift, f0 scale)
        placement of a modelled curve is an integer gather from C. Brightness is the window
        interpolated in time onto tq and along log-frequency onto the axis; in mode="ridge"
        it is the ridge kernel evaluated at each axis frequency."""
        if n in self._comb:
            return self._comb[n]
        tq = np.linspace(self.fit.t_start, self.fit.t_end, n)
        f = self.sp.f
        lf0 = float(np.log(f[0])); L = int(np.log(f[-1] / f[0]) / DLOG) + 1
        lf = lf0 + DLOG * np.arange(L)
        C = np.zeros((n, L), np.float32)
        if self.mode == "brightness":
            nf = len(f)
            Sq = self.sp.sample(np.repeat(tq, nf), np.tile(f, n)).reshape(n, nf)          # time half of the bilinear sample
            Slog = np.stack([np.interp(lf, np.log(f), Sq[k]) for k in range(n)]).astype(np.float32)
            for h in range(1, self.harmonics + 1):
                o = int(round(np.log(h) / DLOG))
                if o < L:
                    C[:, :L - o] += Slog[:, o:]
        else:
            ts, fs, ws = self.ridge_data
            fa = np.exp(lf)
            i = np.rint((tq - ts[0]) / (ts[1] - ts[0])).astype(int)
            inside = (i >= 0) & (i < len(ts))
            ic = np.clip(i, 0, len(ts) - 1)
            for h in range(1, self.harmonics + 1):
                fr = np.interp(tq, ts, fs[h - 1] / h, left=np.nan, right=np.nan)          # the ridge, as a fundamental
                wr = np.where(inside, ws[h - 1][ic], 0.0)
                with np.errstate(invalid="ignore"):
                    z = (fa[None, :] - fr[:, None]) / (self.ridge_sigma * fr[:, None])
                    k = wr[:, None] * np.exp(-0.5 * z * z)
                C += np.nan_to_num(k, nan=0.0).astype(np.float32)
        # zero margins of _PAD columns each side, rows laid end to end: a placement is then a
        # flat index (row * width + column), and anything off the axis reads a zero
        Cp = np.zeros((n, L + 4 * _PAD), np.float32)
        Cp[:, 2 * _PAD:2 * _PAD + L] = C
        self._comb[n] = (tq, Cp.ravel(), lf0, L)
        return self._comb[n]

    def curve_score(self, t_audio, ratio, f0_scale=1.0) -> float:
        """Score of one placement; same convention as sweep_scores."""
        return float(self.sweep_scores(t_audio, ratio, [0.0], [f0_scale])[0, 0])

    def sweep_scores(self, t_audio, ratio, shifts, f0_scales, n: int = 160) -> np.ndarray:
        """Scores for every (shift, f0_scale) combination, (S, K).

        The modelled curve (frequency ratio f/f0 at audio times t_audio) is resampled onto n
        query times of the fit window for each shift, placed on the log-frequency axis of
        comb(n) at log(f0 * ratio), offset by log(f0_scale), and the evidence along it is
        summed and divided by n: the mean over the *whole* fit window, with the parts of the
        window the modelled curve does not cover counting as zero. A flight whose trajectory
        only explains half the measured sweep (it took off, landed or left the data mid-pass)
        earns half the score instead of being judged only where it happens to fit. Fewer than
        10 covered query times scores zero. The result is divided by the number of harmonics.

        Placements are rounded to the axis step DLOG (0.1 % of frequency); against bilinear
        sampling on the linear axis that changes scores by a few hundredths of a dB."""
        tq, Cf, lf0, L = self.comb(n)
        shifts = np.asarray(shifts, float); f0_scales = np.asarray(f0_scales, float)
        ok = np.isfinite(ratio)
        rat = np.interp(tq[None, :] - shifts[:, None], np.asarray(t_audio)[ok], np.asarray(ratio)[ok],
                        left=np.nan, right=np.nan)                                   # (S, n)
        valid = np.isfinite(rat).sum(1) >= 10
        off = np.rint(np.log(f0_scales) / DLOG).astype(np.int64)                   # (K,)
        if np.abs(off).max(initial=0) >= _PAD:
            raise ValueError(f"f0 scales beyond {np.exp(_PAD * DLOG):.3f} are not supported")
        with np.errstate(invalid="ignore", divide="ignore"):
            base = np.rint((np.log(self.fit.f0 * rat) - lf0) / DLOG)               # (S, n) axis column
        base = np.clip(np.nan_to_num(base, nan=-_PAD, posinf=L + _PAD - 1, neginf=-_PAD), -_PAD, L + _PAD - 1)
        flat = base.astype(np.int64) + 2 * _PAD + np.arange(n)[None, :] * (L + 4 * _PAD)
        total = np.take(Cf, flat[:, None, :] + off[None, :, None]).sum(axis=2, dtype=np.float64) / n   # (S, K)
        total[~valid] = 0.0
        return total / self.harmonics

    def emission_span(self, t0_map: float, max_shift: float = 16.0) -> tuple[float, float]:
        """Unix emission times a modelled curve needs to cover this window under the clock
        audio = unix + t0_map, for shifts up to max_shift: the fit window, the shifts, and up
        to 80 s of propagation delay (27 km)."""
        return self.fit.t_start - t0_map - max_shift - 80.0, self.fit.t_end - t0_map + max_shift + 1.0


def modelled_curve(track: Track, obs_enu, c: float):
    """(t_obs unix, ratio f/f0) the observer at obs_enu hears from a sampled track: the range
    rate at each emission time from the track's own velocity, heard r / c later."""
    d = track.pos - np.asarray(obs_enu, float)
    r = np.linalg.norm(d, axis=1)
    rr = (d * track.vel).sum(1) / np.maximum(r, 1e-6)
    return track.t + r / c, c / (c + rr)


def locked_score(ev1: EventWindow, ev2: EventWindow, tr1: Track | None, tr2: Track | None, obs_enu,
                 t_heard1: float, c: float | None = None, shifts=SHIFTS, f0_scales=F0_SCALES):
    """Two-event locked scoring.

    The recording's clock is set so that event 1's modelled heard CPA (t_heard1, unix) lands on
    the measured one. Event 1's placement is then swept over small time shifts and f0 scales
    (both are measurements). Event 2 is placed with the *same* time shift, so its timing
    relative to event 1 is fixed purely by flight geometry. A track that is None (the flight
    does not cover the window) scores zero. Returns (total, s1, s2, best_shift)."""
    c = speed_of_sound() if c is None else c
    t0_map = ev1.fit.t_cpa - t_heard1
    if tr1 is None:
        return 0.0, 0.0, 0.0, 0.0
    to1, r1 = modelled_curve(tr1, obs_enu, c)
    sc1 = ev1.sweep_scores(to1 + t0_map, r1, shifts, f0_scales)
    i, j = np.unravel_index(np.argmax(sc1), sc1.shape)
    s1, sh = float(sc1[i, j]), float(shifts[i])
    if tr2 is None:
        return s1, s1, 0.0, sh
    to2, r2 = modelled_curve(tr2, obs_enu, c)
    s2 = float(ev2.sweep_scores(to2 + t0_map, r2, [sh], f0_scales).max())
    return s1 + s2, s1, s2, sh


def single_score(ev: EventWindow, tr: Track | None, obs_enu, t_heard: float, c: float | None = None,
                 shifts=SHIFTS, f0_scales=F0_SCALES) -> tuple[float, float]:
    """One event on the clock of its own heard CPA: (best score, its shift)."""
    c = speed_of_sound() if c is None else c
    if tr is None:
        return 0.0, 0.0
    to, r = modelled_curve(tr, obs_enu, c)
    sc = ev.sweep_scores(to + (ev.fit.t_cpa - t_heard), r, shifts, f0_scales)
    i, j = np.unravel_index(np.argmax(sc), sc.shape)
    return float(sc[i, j]), float(shifts[i])

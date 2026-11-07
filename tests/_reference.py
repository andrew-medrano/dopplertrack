"""Reference implementations the fast code is checked against (not used by the package).

linear_sweep_scores is the dopplertrack 0.3 scorer: the modelled curve placed on the linear
frequency axis and the window sampled bilinearly, every harmonic separately."""
import numpy as np


def linear_sweep_scores(w, t_audio, ratio, shifts, f0_scales, n=160):
    """EventWindow.sweep_scores as dopplertrack 0.3 computed it (brightness or ridge mode)."""
    tq = np.linspace(w.fit.t_start, w.fit.t_end, n)
    shifts = np.asarray(shifts, float); f0_scales = np.asarray(f0_scales, float)
    ok = np.isfinite(ratio)
    rat = np.interp(tq[None, :] - shifts[:, None], t_audio[ok], ratio[ok], left=np.nan, right=np.nan)   # (S, n)
    valid = np.isfinite(rat).sum(1) >= 10
    total = np.zeros((len(shifts), len(f0_scales)))
    T = np.broadcast_to(tq[None, None, :], (len(shifts), len(f0_scales), n))
    for h in range(1, w.harmonics + 1):
        F = h * w.fit.f0 * f0_scales[None, :, None] * rat[:, None, :]           # (S, K, n)
        if w.mode == "brightness":
            v = w.sp.sample(T.ravel(), np.nan_to_num(F.ravel(), nan=-1.0)).reshape(F.shape)
        else:
            ts, fs, ws = w.ridge_data
            i = np.rint((T - ts[0]) / (ts[1] - ts[0])).astype(int)
            inside = (i >= 0) & (i < len(ts))
            fr = np.interp(T, ts, fs[h - 1], left=np.nan, right=np.nan)
            wr = np.where(inside, ws[h - 1][np.clip(i, 0, len(ts) - 1)], 0.0)
            with np.errstate(invalid="ignore"):
                z = (F - fr) / (w.ridge_sigma * fr)
                v = np.nan_to_num(wr * np.exp(-0.5 * z * z), nan=0.0)
        total += np.where(np.isfinite(F), v, 0.0).sum(axis=2) / n
    total[~valid] = 0.0
    return total / w.harmonics

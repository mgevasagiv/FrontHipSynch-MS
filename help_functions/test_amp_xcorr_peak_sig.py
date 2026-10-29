#!/usr/bin/env python3
"""Smoke test mirroring help_functions/amp_xcorr_peak_sig.m (no MATLAB required)."""

from __future__ import annotations

import numpy as np


def xcorr_coeff(a: np.ndarray, b: np.ndarray, max_lag: int) -> tuple[np.ndarray, np.ndarray]:
    """MATLAB-compatible xcorr(a,b): positive lag => b delayed relative to a."""
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    # np.correlate uses the opposite lag sign vs MATLAB xcorr; flip via correlate(b, a)
    full = np.correlate(b, a, mode="full")
    mid = len(a) - 1
    xc = full[mid - max_lag : mid + max_lag + 1]
    denom = np.sqrt(np.sum(a**2) * np.sum(b**2))
    if denom == 0:
        xc = np.zeros_like(xc)
    else:
        xc = xc / denom
    lags = np.arange(-max_lag, max_lag + 1)
    return xc, lags


def amp_xcorr_peak_sig(
    a: np.ndarray,
    b: np.ndarray,
    *,
    fs: float = 512.0,
    max_lag_sec: float = 0.2,
    n_perm: int = 1000,
    alpha: float = 0.05,
    seed: int | None = None,
) -> dict:
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    assert a.shape == b.shape
    rng = np.random.default_rng(seed)

    a0 = a - np.nanmean(a)
    b0 = b - np.nanmean(b)
    a0 = np.nan_to_num(a0, nan=0.0)
    b0 = np.nan_to_num(b0, nan=0.0)

    max_lag = min(int(round(max_lag_sec * fs)), len(a0) - 1)
    xc, lags = xcorr_coeff(a0, b0, max_lag)
    peak_idx = int(np.argmax(xc))
    peak_value = float(xc[peak_idx])
    peak_lag = int(lags[peak_idx])

    null_peaks = np.empty(n_perm)
    for i in range(n_perm):
        bp = rng.permutation(b0)
        xc_p, _ = xcorr_coeff(a0, bp, max_lag)
        null_peaks[i] = float(np.max(xc_p))

    null_thresh = float(np.percentile(null_peaks, 100 * (1 - alpha)))
    p_peak = (1 + np.sum(null_peaks >= peak_value)) / (1 + n_perm)
    return {
        "peakValue": peak_value,
        "peakLagSec": peak_lag / fs,
        "nullThresh": null_thresh,
        "peakSignificant": peak_value > null_thresh,
        "p_peak": float(p_peak),
    }


def main() -> None:
    fs = 512.0
    t = np.arange(0, 1.5, 1 / fs)
    rng = np.random.default_rng(0)
    a = np.sin(2 * np.pi * 3 * t) + 0.1 * rng.normal(size=t.shape)
    true_lag_sec = 0.08
    shift = int(round(true_lag_sec * fs))
    b = np.concatenate([np.zeros(shift), a[:-shift]]) + 0.1 * rng.normal(size=t.shape)

    out = amp_xcorr_peak_sig(a, b, fs=fs, max_lag_sec=0.2, n_perm=500, seed=42)
    assert out["peakSignificant"], out
    assert abs(out["peakLagSec"] - true_lag_sec) <= 2 / fs, out

    out_null = amp_xcorr_peak_sig(
        rng.normal(size=t.shape),
        rng.normal(size=t.shape),
        fs=fs,
        max_lag_sec=0.2,
        n_perm=500,
        seed=7,
    )
    print(
        f"OK lagged: lag={out['peakLagSec']:+.4f}s p={out['p_peak']:.4g} | "
        f"null sig={out_null['peakSignificant']} p={out_null['p_peak']:.4g}"
    )


if __name__ == "__main__":
    main()

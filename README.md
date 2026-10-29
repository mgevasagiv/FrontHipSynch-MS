# FrontHipSynchron-MS
Code to reproduce 'Fronto-Hippocampal Synchronization' manuscript analysis and stats

## Amplitude cross-correlation peak significance (Figure 2D)

- `help_functions/amp_xcorr_peak_sig.m` — Adhikari-style mean-subtracted `xcorr` over ±0.2 s, with a 1000-permutation null of peak values (peak significant if > 95th percentile of null).
- `MAZE10_plotTFR_OFC_axis.m` — loads OFC block-1/block-2 population TFRs from `MAZE10_calc_TFR_event`, builds mid-gamma time courses, runs the shuffle test for events `X` / `G` / `N`. Set `doDemo = true` for a synthetic self-check without data.
- `help_functions/test_amp_xcorr_peak_sig.py` — Python smoke test of the same null logic.

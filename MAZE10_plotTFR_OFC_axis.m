% MAZE10_plotTFR_OFC_axis
% Block-1 vs block-2 OFC mid-gamma amplitude cross-correlation (Figure 2D)
%
% Loads population TFRs written by MAZE10_calc_TFR_event('OFC') for improved
% routes (cond 0 = block 1, cond 1 = block 2), builds peri-event mid-gamma
% time courses (freq-average over params.mgamma_range, electrode-average),
% then ascertains cross-correlation peak significance with the manuscript
% shuffle procedure (help_functions/amp_xcorr_peak_sig.m):
%   - mean-subtract envelopes
%   - xcorr over lags +/- 0.2 s
%   - 1000 random permutations -> null of peak values
%   - peak significant if observed peak > 95th percentile of null
%
% Also reports lag-0 corrcoef r/p as in the manuscript text (descriptive
% only; peak inference uses the shuffle test).
%
% Prerequisites:
%   CenteredAvgEvents_allElectrodesSubj_OFC_0.mat  (block 1, improved)
%   CenteredAvgEvents_allElectrodesSubj_OFC_1.mat  (block 2, improved)
% under dirs.TF/TFR_grp/OFC/population/ (see maze_set_path).
%
% Optional: set doDemo = true to run a synthetic self-check without data.

clear all; close all;

doDemo = false;          % true -> synthetic signals only (no MAT load)
nPerm = 1000;
maxLagSec = 0.2;
alpha = 0.05;
srate = 512;             % matches MAZE10_calc_TFR_event downsampled rate
periSec = 0.75;          % params.timeBeforeAfterEventRipSpec in MAZE10
events_of_interest = {'X', 'G', 'N'};  % maze onset, goal-reveal, open path

addpath(genpath(fileparts(mfilename('fullpath'))));
if exist('maze_set_path', 'file') == 2 && ~doDemo
    try
        maze_set_path;
    catch
        warning('maze_set_path failed; set savedir manually if needed.');
    end
end

%% Synthetic demo (no data required)
if doDemo
    demo_amp_xcorr_peak_sig();
    return
end

%% Paths / load block-wise population TFRs
if exist('dirs', 'var') && isfield(dirs, 'TF')
    popdir = fullfile(dirs.TF, 'TFR_grp', 'OFC', 'population');
    plotdir = fullfile(dirs.TF, 'TFR_grp', 'population');
else
    popdir = fullfile('E:', 'MAZE_1.0', 'analysis', 'TF', 'TFR_grp', 'OFC', 'population');
    plotdir = fullfile('E:', 'MAZE_1.0', 'analysis', 'TF', 'TFR_grp', 'population');
end
if ~exist(plotdir, 'dir'); mkdir(plotdir); end

f_b1 = fullfile(popdir, 'CenteredAvgEvents_allElectrodesSubj_OFC_0.mat');
f_b2 = fullfile(popdir, 'CenteredAvgEvents_allElectrodesSubj_OFC_1.mat');
assert(exist(f_b1, 'file') == 2, 'Missing block-1 file: %s', f_b1);
assert(exist(f_b2, 'file') == 2, 'Missing block-2 file: %s', f_b2);

B1 = load(f_b1, 'meanTFREv', 'params', 'evlist', 'subjects');
B2 = load(f_b2, 'meanTFREv', 'params', 'evlist', 'subjects');

if isfield(B1.params, 'mgamma_range')
    gamma_range = B1.params.mgamma_range;
else
    gamma_range = [50 100];
end
if isfield(B1.params, 'freqRangeForAvgSpec')
    freqs = B1.params.freqRangeForAvgSpec;
else
    freqs = 1:250;
end
fMask = freqs >= gamma_range(1) & freqs <= gamma_range(2);
tvec = -periSec:1/srate:periSec;

%% Per-event cross-correlation + shuffle peak test
results = struct('event', {}, 'nTime', {}, 'r0', {}, 'p0', {}, ...
    'peakValue', {}, 'peakLagSec', {}, 'nullThresh', {}, ...
    'peakSignificant', {}, 'p_peak', {}, 'nPerm', {});

fprintf('\n=== OFC mid-gamma block1-vs-block2 amplitude xcorr (nPerm=%d) ===\n', nPerm);
for iE = 1:numel(events_of_interest)
    ev = events_of_interest{iE};
    if ~isfield(B1.meanTFREv, ev) || ~isfield(B2.meanTFREv, ev)
        warning('Event %s missing from meanTFREv; skipping.', ev);
        continue
    end

    tc1 = gamma_timecourse(B1.meanTFREv.(ev), fMask, tvec);
    tc2 = gamma_timecourse(B2.meanTFREv.(ev), fMask, tvec);
    n = min(numel(tc1), numel(tc2));
    tc1 = tc1(1:n);
    tc2 = tc2(1:n);
    tUse = tvec(1:n);

    out = amp_xcorr_peak_sig(tc1, tc2, ...
        'fs', srate, 'maxLagSec', maxLagSec, 'nPerm', nPerm, ...
        'alpha', alpha, 'seed', 1);

    results(iE).event = ev; %#ok<*SAGROW>
    results(iE).nTime = n;
    results(iE).r0 = out.r0;
    results(iE).p0 = out.p0;
    results(iE).peakValue = out.peakValue;
    results(iE).peakLagSec = out.peakLagSec;
    results(iE).nullThresh = out.nullThresh;
    results(iE).peakSignificant = out.peakSignificant;
    results(iE).p_peak = out.p_peak;
    results(iE).nPerm = out.nPerm;

    fprintf(['%s: peak=%.3f at lag=%+.3f s | shuffle p=%.4g (sig@%.0f%%=%d) | ' ...
        'lag0 r=%.3f p_corrcoef=%.3g\n'], ...
        ev, out.peakValue, out.peakLagSec, out.p_peak, 100*(1-alpha), ...
        out.peakSignificant, out.r0, out.p0);

    % --- figure panel ---
    fig_name = sprintf('OFC_gamma_xcorr_blk12_%s', ev);
    figure('Name', fig_name, 'Color', 'w');

    subplot(2, 2, 1); hold on
    plot(tUse, tc1, 'k', 'LineWidth', 1.5);
    plot(tUse, tc2, 'Color', [0.7 0.2 0.2], 'LineWidth', 1.5);
    xline(0, 'k:');
    xlabel('Time from event (s)'); ylabel('Mid-gamma power (a.u.)');
    legend({'Block 1', 'Block 2'}, 'Location', 'best');
    title(sprintf('%s gamma time courses', ev));

    subplot(2, 2, 2); hold on
    plot(out.lagsSec, out.xcorr, 'k', 'LineWidth', 1.2);
    yline(out.nullThresh, 'r--', 'LineWidth', 1);
    plot(out.peakLagSec, out.peakValue, 'ro', 'MarkerFaceColor', 'r');
    xlabel('Lag (s)'); ylabel('xcorr (coeff)');
    title(sprintf('xcorr peak @ %+.3f s', out.peakLagSec));
    legend({'xcorr', sprintf('null %.0f%%', 100*(1-alpha)), 'peak'}, ...
        'Location', 'best');

    subplot(2, 2, 3); hold on
    histogram(out.nullPeaks, 40, 'Normalization', 'pdf', 'FaceColor', [0.7 0.7 0.7]);
    xline(out.peakValue, 'r', 'LineWidth', 2);
    xline(out.nullThresh, 'r--', 'LineWidth', 1);
    xlabel('Peak xcorr'); ylabel('Density');
    title(sprintf('Shuffle null (p_{peak}=%.3g)', out.p_peak));

    subplot(2, 2, 4); axis off
    txt = {
        sprintf('Event: %s', ev)
        sprintf('Peak value: %.4f', out.peakValue)
        sprintf('Peak lag: %+.4f s', out.peakLagSec)
        sprintf('Null %.0f%% thresh: %.4f', 100*(1-alpha), out.nullThresh)
        sprintf('Peak significant: %d', out.peakSignificant)
        sprintf('Shuffle p_peak: %.4g', out.p_peak)
        sprintf('corrcoef lag0 r: %.4f', out.r0)
        sprintf('corrcoef lag0 p: %.4g  (descriptive only)', out.p0)
        sprintf('nPerm: %d', out.nPerm)
        };
    text(0.05, 0.5, txt, 'FontName', 'FixedWidth', 'FontSize', 11, ...
        'VerticalAlignment', 'middle');

    saveas(gcf, fullfile(plotdir, fig_name), 'png');
end

save(fullfile(plotdir, 'OFC_gamma_xcorr_blk12_results.mat'), ...
    'results', 'nPerm', 'maxLagSec', 'alpha', 'gamma_range', 'srate');

fprintf('\nSaved results to %s\n', fullfile(plotdir, 'OFC_gamma_xcorr_blk12_results.mat'));


%% -------- helpers --------
function tc = gamma_timecourse(meanTFR, fMask, tvec)
% meanTFR: freq x time x electrodes (from MAZE10_calc_TFR_event)
if ndims(meanTFR) < 2
    error('meanTFR must be at least freq x time');
end
nFreq = size(meanTFR, 1);
nTime = size(meanTFR, 2);
if numel(fMask) ~= nFreq
    error('Frequency mask length (%d) != nFreq (%d)', numel(fMask), nFreq);
end
band = meanTFR(fMask, :, :);                 % nBand x time x elec
tc = squeeze(nanmean(nanmean(band, 1), 3));  % time x 1
tc = tc(:);
if numel(tc) ~= numel(tvec)
    % Allow slight length mismatch from inclusive endpoint sampling
    n = min(numel(tc), numel(tvec));
    tc = tc(1:n);
end
end


function demo_amp_xcorr_peak_sig()
% Synthetic self-check: B is a lagged copy of A -> significant positive lag peak
fs = 512;
t = (0:1/fs:1.5)';
A = sin(2*pi*3*t) + 0.1*randn(size(t));
trueLagSec = 0.08;
shift = round(trueLagSec * fs);
B = [zeros(shift, 1); A(1:end-shift)] + 0.1*randn(size(t));

out = amp_xcorr_peak_sig(A, B, 'fs', fs, 'maxLagSec', 0.2, ...
    'nPerm', 500, 'seed', 42);

fprintf('DEMO: true lag = +%.3f s | recovered lag = %+.3f s | peakSignificant = %d | p_peak = %.4g\n', ...
    trueLagSec, out.peakLagSec, out.peakSignificant, out.p_peak);
assert(out.peakSignificant == 1, 'Demo failed: expected significant peak');
assert(abs(out.peakLagSec - trueLagSec) <= 2/fs, 'Demo failed: lag recovery');

% Unrelated noise should usually fail significance
outNull = amp_xcorr_peak_sig(randn(size(A)), randn(size(A)), ...
    'fs', fs, 'maxLagSec', 0.2, 'nPerm', 500, 'seed', 7);
fprintf('DEMO null: peakSignificant = %d | p_peak = %.4g\n', ...
    outNull.peakSignificant, outNull.p_peak);
fprintf('Demo OK.\n');
end

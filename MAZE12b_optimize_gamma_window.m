% MAZE12b_optimize_gamma_window
% -----------------------------------------------------------------------------
% Data-driven selection of the post-stimulus OFC gamma (50-100 Hz) analysis
% window, responding to the reviewer concern that Methods said [0.15 0.35]
% while code/figure used [0.2 0.35].
%
% Pipeline
%   1) Load improved-trial OFC TFRs (block1=cond0, block2=cond1) from MAZE10
%   2) Build gamma time courses γ(t) for events X, N, G, preG
%   3) Cluster-based permutation of (high-load − open-path) contrast over time
%   4) Define optimal_window from the earliest significant post-stimulus cluster
%   5) Confirmatory GCP LMEs in that window (incl. single block-2 model)
%   6) Robustness grid across nearby windows
%
% Inputs (written by MAZE10_calc_TFR_event('OFC')):
%   E:\MAZE_1.0\analysis\TF\TFR_grp\OFC\population\
%     CenteredAvgEvents_allElectrodesSubj_OFC_0.mat           % block1 improved
%     CenteredAvgEvents_allElectrodesSubj_OFC_1.mat           % block2 improved
%     CenteredAvgEvents_allElectrodesSubj_OFC_*_TFRMAT.mat    % optional trial TFRs
%
% Outputs:
%   ...\population\gamma_window_optimization\
%     optimal_gamma_window.mat / .csv
%     cluster_timecourse.png
%     gcp_confirmatory_stats.mat
%     robustness_grid.csv
%     fig2CF_unified_block2.png
%
% Run:  MAZE12b_optimize_gamma_window

clear; close all; clc;

%% ---------------- paths / params ----------------
try
    maze_set_path;
catch
    warning('maze_set_path failed; using hardcoded MAZE_1.0 roots');
end
if ~exist('dirs', 'var') || ~isfield(dirs, 'TF')
    dirs.TF = 'E:\MAZE_1.0\analysis\TF\';
end

popdir  = fullfile(dirs.TF, 'TFR_grp', 'OFC', 'population');
outdir  = fullfile(popdir, 'gamma_window_optimization');
if ~exist(outdir, 'dir'); mkdir(outdir); end

events_hl = {'X','G','preG'};   % high computational load
event_bl  = 'N';                % open-path baseline
events_all = [{'N'}, events_hl];

fb_mgamma = [50 100];
legacy_window = [0.2 0.35];     % what the code actually used
text_window   = [0.15 0.35];    % what the manuscript text said (typo)
search_lim    = [0 0.60];       % only search post-stimulus clusters here
nperm         = 2000;
cluster_alpha = 0.05;           % pointwise threshold for cluster formation
cluster_p     = 0.05;           % cluster-level significance
rng(42);

fprintf('Population dir: %s\n', popdir);

%% ---------------- load block 1 / 2 population mats ----------------
blk_conds = [0 1];   % improved routes: block1, block2
blk_ids   = [1 2];
S = struct([]);

for ib = 1:numel(blk_conds)
    cond = blk_conds(ib);
    fpop = fullfile(popdir, sprintf('CenteredAvgEvents_allElectrodesSubj_OFC_%d.mat', cond));
    assert(isfile(fpop), 'Missing population file: %s', fpop);
    mm = matfile(fpop);
    S(ib).cond = cond;
    S(ib).blk  = blk_ids(ib);
    S(ib).params = mm.params;
    S(ib).subjects = mm.subjects;
    S(ib).evlist = mm.evlist;
    S(ib).meanTFREv = mm.meanTFREv;
    S(ib).X_subj = mm.X_subj;
    S(ib).X_elec = mm.X_elec;
    S(ib).nEvents = mm.nEvents;
    % optional trial-resolved TFRs
    ftfr = fullfile(popdir, sprintf('CenteredAvgEvents_allElectrodesSubj_OFC_%d_TFRMAT.mat', cond));
    if isfile(ftfr)
        mt = matfile(ftfr);
        S(ib).ALL_TFR_MAT = mt.ALL_TFR_MAT;
        fprintf('Loaded trial TFRs for cond %d\n', cond);
    else
        S(ib).ALL_TFR_MAT = [];
        fprintf('No TFRMAT for cond %d — electrode-level discovery only\n', cond);
    end
end

params = S(1).params;
srate  = 512;
if isfield(params, 'timeBeforeAfterEventRipSpec')
    twin = params.timeBeforeAfterEventRipSpec;
else
    twin = 0.75;
end
tvec = -twin:1/srate:twin;
freq = params.freqRangeForAvgSpec(:);
fidx = find(freq >= fb_mgamma(1) & freq <= fb_mgamma(2));
assert(~isempty(fidx), 'No frequencies in gamma band');

%% ---------------- build subject × event × block gamma time courses -------
% Discovery unit = subject (mean over that subject's OFC electrodes)
subjects = S(1).subjects;
nSub = numel(subjects);
nT = numel(tvec);
% Gtc(subj, eventIndex, blkIndex, time)
Gtc = nan(nSub, numel(events_all), numel(S), nT);

for ib = 1:numel(S)
    for ie = 1:numel(events_all)
        ev = events_all{ie};
        assert(isfield(S(ib).meanTFREv, ev), 'Missing meanTFREv.%s', ev);
        M = S(ib).meanTFREv.(ev);          % freq × time × elec
        subj_e = S(ib).X_subj.(ev)(:);     % elec × 1 (subject index into subjects list)
        assert(size(M,3) == numel(subj_e), 'meanTFREv / X_subj size mismatch for %s', ev);
        g_elec = squeeze(nanmean(M(fidx, :, :), 1));   % time × elec
        if isvector(g_elec); g_elec = g_elec(:); end
        for isub = 1:nSub
            hit = (subj_e == isub);
            if ~any(hit); continue; end
            Gtc(isub, ie, ib, :) = nanmean(g_elec(:, hit), 2);
        end
    end
end

% High-load = mean of X,G,preG ; Open-path = N
idxN = find(strcmp(events_all, 'N'));
idxHL = find(ismember(events_all, events_hl));
% Contrast per subject × block × time
C_sb = nan(nSub, numel(S), nT);
for ib = 1:numel(S)
    gN  = squeeze(Gtc(:, idxN, ib, :));                 % subj × time
    gHL = squeeze(nanmean(Gtc(:, idxHL, ib, :), 2));    % subj × time
    C_sb(:, ib, :) = gHL - gN;
end
% Pool blocks for discovery (primary): subject mean contrast across blocks
C_subj = squeeze(nanmean(C_sb, 2));   % subj × time

%% ---------------- observed t-map + clusters ----------------
[t_obs, clusters_obs] = contrast_timecourse_stats(C_subj, tvec, search_lim, cluster_alpha);

%% ---------------- cluster permutation ----------------
max_cluster_mass_perm = zeros(nperm, 1);
nSubEff = sum(any(isfinite(C_subj), 2));
fprintf('Cluster permutation (%d perms), nSub with data = %d ...\n', nperm, nSubEff);
for ip = 1:nperm
    flips = sign(randn(nSub, 1));
    flips(flips == 0) = 1;
    C_perm = C_subj .* flips;
    [~, clusters_p] = contrast_timecourse_stats(C_perm, tvec, search_lim, cluster_alpha);
    if isempty(clusters_p)
        max_cluster_mass_perm(ip) = 0;
    else
        max_cluster_mass_perm(ip) = max(abs([clusters_p.mass]));
    end
    if mod(ip, 500) == 0
        fprintf('  perm %d/%d\n', ip, nperm);
    end
end
mass_thresh = prctile(max_cluster_mass_perm, 100*(1-cluster_p));

sig_clusters = clusters_obs;
for ic = 1:numel(sig_clusters)
    sig_clusters(ic).p_corr = mean(max_cluster_mass_perm >= abs(sig_clusters(ic).mass));
    sig_clusters(ic).significant = sig_clusters(ic).p_corr < cluster_p;
end

% Choose earliest significant POSITIVE post-stimulus cluster as optimal window
pos_sig = sig_clusters([sig_clusters.significant] & [sig_clusters.mean_t] > 0);
if isempty(pos_sig)
    warning('No significant positive post-stimulus cluster. Falling back to legacy [0.2 0.35].');
    optimal_window = legacy_window;
    chosen = struct('t_onset', legacy_window(1), 't_offset', legacy_window(2), ...
        'p_corr', NaN, 'mass', NaN, 'mean_t', NaN, 'significant', false, 'fallback', true);
else
    [~, ord] = sort([pos_sig.t_onset]);
    chosen = pos_sig(ord(1));
    chosen.fallback = false;
    % Snap to 10 ms grid for a clean Methods number; keep >=50 ms wide
    t1 = floor(chosen.t_onset*100)/100;
    t2 = ceil(chosen.t_offset*100)/100;
    if (t2 - t1) < 0.05
        t2 = t1 + 0.05;
    end
    optimal_window = [t1 t2];
end

fprintf('\n=== Optimal gamma window ===\n');
fprintf('Selected: [%.3f %.3f] s  (legacy code=[%.2f %.2f], text typo=[%.2f %.2f])\n', ...
    optimal_window(1), optimal_window(2), legacy_window(1), legacy_window(2), text_window(1), text_window(2));
if isfield(chosen,'p_corr')
    fprintf('Cluster p_corr=%.4f, mass=%.2f, mean_t=%.2f\n', chosen.p_corr, chosen.mass, chosen.mean_t);
end

%% ---------------- plot discovery time course ----------------
fig = figure('Color','w','Position',[100 100 900 420]);
mu = nanmean(C_subj, 1);
se = nanstd(C_subj, 0, 1) ./ sqrt(sum(isfinite(C_subj),1));
hold on;
% shade significant clusters
yl = [min(mu-se)*1.2, max(mu+se)*1.2];
if isempty(yl) || any(~isfinite(yl)); yl = [-1 1]; end
for ic = 1:numel(sig_clusters)
    if ~sig_clusters(ic).significant; continue; end
    x = [sig_clusters(ic).t_onset, sig_clusters(ic).t_offset, sig_clusters(ic).t_offset, sig_clusters(ic).t_onset];
    y = [yl(1) yl(1) yl(2) yl(2)];
    patch(x, y, [1 0.85 0.4], 'EdgeColor','none', 'FaceAlpha', 0.35);
end
plot(tvec, mu, 'k', 'LineWidth', 1.5);
plot(tvec, mu+se, 'Color',[0.5 0.5 0.5]);
plot(tvec, mu-se, 'Color',[0.5 0.5 0.5]);
plot([0 0], yl, 'k--');
plot(optimal_window, [0 0], 'r-', 'LineWidth', 4);
xlim([-0.2 0.6]); ylim(yl);
xlabel('Time from behavioral event (s)');
ylabel('OFC \gamma_{50-100} high-load − open-path (subj mean)');
title(sprintf('Cluster-defined window [%.2f %.2f] s (nperm=%d)', optimal_window(1), optimal_window(2), nperm));
saveas(fig, fullfile(outdir, 'cluster_timecourse.png'));
saveas(fig, fullfile(outdir, 'cluster_timecourse.fig'));

%% ---------------- confirmatory GCP tables ----------------
% Electrode-level GCP in optimal / legacy / text windows
windows = struct( ...
    'optimal', optimal_window, ...
    'legacy_code', legacy_window, ...
    'text_typo', text_window, ...
    'w_015_040', [0.15 0.40], ...
    'w_020_040', [0.20 0.40], ...
    'w_015_030', [0.15 0.30]);

wnames = fieldnames(windows);
robust_rows = {};
T_opt = table();

for iw = 1:numel(wnames)
    w = windows.(wnames{iw});
    T = build_gcp_table_electrode(S, events_all, fidx, tvec, w);
    % Full factorial (both blocks)
    T.event = categorical(T.event, events_all);
    lme = fitlme(T, 'GCP ~ event + blk + event:blk + (1|subj)');
    an = anova(lme);
    % Block-2 unified model (answers Fig 2C/2F reviewer point)
    T2 = T(T.blk == 2, :);
    lme2 = fitlme(T2, 'GCP ~ event + (1|subj)');
    an2 = anova(lme2);
    % Planned contrasts vs N within block 2
    contrast_p = nan(1, numel(events_hl));
    for ie = 1:numel(events_hl)
        Ts = T2(T2.event == event_bl | T2.event == events_hl{ie}, :);
        Ts.event = removecats(Ts.event);
        try
            lm = fitlme(Ts, 'GCP ~ event + (1|subj)');
            contrast_p(ie) = lm.Coefficients.pValue(2);
        catch
            contrast_p(ie) = NaN;
        end
    end
    % Holm correction across the 3 contrasts
    contrast_p_holm = holm_correct(contrast_p);

    robust_rows(end+1, :) = {wnames{iw}, w(1), w(2), ... %#ok<AGROW>
        an.pValue(strcmp(an.Term,'event')), ...
        an.pValue(strcmp(an.Term,'blk')), ...
        an.pValue(strcmp(an.Term,'event:blk')), ...
        an2.pValue(strcmp(an2.Term,'event')), ...
        contrast_p(1), contrast_p(2), contrast_p(3), ...
        contrast_p_holm(1), contrast_p_holm(2), contrast_p_holm(3)};

    if strcmp(wnames{iw}, 'optimal')
        T_opt = T;
        stats_opt.lme_full = lme;
        stats_opt.anova_full = an;
        stats_opt.lme_block2 = lme2;
        stats_opt.anova_block2 = an2;
        stats_opt.contrast_events = events_hl;
        stats_opt.contrast_p = contrast_p;
        stats_opt.contrast_p_holm = contrast_p_holm;
        fprintf('\n=== Confirmatory LME (optimal window) ===\n');
        disp(an);
        fprintf('Block-2 unified event effect p = %.3g\n', an2.pValue(strcmp(an2.Term,'event')));
        for ie = 1:numel(events_hl)
            fprintf('  Block2 %s vs N: p=%.3g (Holm=%.3g)\n', events_hl{ie}, contrast_p(ie), contrast_p_holm(ie));
        end
    end
end

robustTbl = cell2table(robust_rows, 'VariableNames', ...
    {'window','t1','t2','p_event','p_blk','p_event_blk','p_event_block2', ...
     'p_X_vs_N','p_G_vs_N','p_preG_vs_N','p_X_vs_N_holm','p_G_vs_N_holm','p_preG_vs_N_holm'});
disp(robustTbl);
writetable(robustTbl, fullfile(outdir, 'robustness_grid.csv'));

%% ---------------- Fig 2C/F style: one shared open-path, unified block-2 ----
fig2 = figure('Color','w','Position',[100 100 900 420]);
T2 = T_opt(T_opt.blk == 2, :);
% Panel-like distributions for X,N,G and N,preG,G from THE SAME open-path rows
plot_events_C = {'X','N','G'};
plot_events_F = {'N','preG','G'};
subplot(1,2,1); plot_event_violins(T2, plot_events_C); title('Block 2: outset / open-path / goal (unified N)');
subplot(1,2,2); plot_event_violins(T2, plot_events_F); title('Block 2: open-path / pre-goal / goal (same N)');
saveas(fig2, fullfile(outdir, 'fig2CF_unified_block2.png'));
saveas(fig2, fullfile(outdir, 'fig2CF_unified_block2.fig'));

%% ---------------- optional trial-level GCP if TFRMAT present -------------
trial_stats = struct();
if ~isempty(S(1).ALL_TFR_MAT) && ~isempty(S(2).ALL_TFR_MAT)
    try
        Tt = build_gcp_table_trial(S, events_all, fidx, tvec, optimal_window);
        Tt.event = categorical(Tt.event, events_all);
        trial_stats.lme_full = fitlme(Tt, 'GCP ~ event + blk + event:blk + (1|subj)');
        trial_stats.anova_full = anova(trial_stats.lme_full);
        Tt2 = Tt(Tt.blk == 2, :);
        trial_stats.lme_block2 = fitlme(Tt2, 'GCP ~ event + (1|subj)');
        trial_stats.anova_block2 = anova(trial_stats.lme_block2);
        fprintf('\n=== Trial-level confirmatory (ALL_TFR_MAT) ===\n');
        disp(trial_stats.anova_full);
        writetable(Tt, fullfile(outdir, 'gcp_trial_optimal_window.csv'));
    catch ME
        warning('Trial-level confirmatory failed: %s', ME.message);
    end
end

%% ---------------- save ----------------
save(fullfile(outdir, 'optimal_gamma_window.mat'), ...
    'optimal_window', 'legacy_window', 'text_window', 'chosen', 'sig_clusters', ...
    't_obs', 'tvec', 'C_subj', 'mass_thresh', 'nperm', 'cluster_alpha', 'cluster_p', ...
    'stats_opt', 'robustTbl', 'T_opt', 'trial_stats', 'fb_mgamma', 'search_lim', '-v7.3');

opt_csv = table(optimal_window(1), optimal_window(2), legacy_window(1), legacy_window(2), ...
    text_window(1), text_window(2), chosen.fallback, ...
    'VariableNames', {'opt_t1','opt_t2','legacy_t1','legacy_t2','text_t1','text_t2','used_fallback'});
writetable(opt_csv, fullfile(outdir, 'optimal_gamma_window.csv'));

fprintf('\nWrote results to %s\n', outdir);
fprintf('Done.\n');

%% ========================= local functions =========================
function [tvec_out, clusters] = contrast_timecourse_stats(C_subj, tvec, search_lim, alpha)
% C_subj: subj × time contrast (high-load − N)
    nSub = size(C_subj, 1);
    nT = size(C_subj, 2);
    tstat = nan(1, nT);
    for it = 1:nT
        x = C_subj(:, it);
        x = x(isfinite(x));
        if numel(x) < 3; continue; end
        [~, ~, ~, st] = ttest(x);  % H0: mean==0
        tstat(it) = st.tstat;
    end
    % two-tailed critical t
    df = max(nSub - 1, 1);
    tcrit = tinv(1 - alpha/2, df);
    in_search = tvec >= search_lim(1) & tvec <= search_lim(2);
    sig = (abs(tstat) >= tcrit) & in_search;
    clusters = form_clusters(sig, tstat, tvec);
    tvec_out = tstat;
end

function clusters = form_clusters(sig, tstat, tvec)
    clusters = struct('idx', {}, 't_onset', {}, 't_offset', {}, 'mass', {}, 'mean_t', {});
    if ~any(sig); return; end
    d = diff([false, sig, false]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    for ic = 1:numel(starts)
        idx = starts(ic):ends(ic);
        clusters(ic).idx = idx;
        clusters(ic).t_onset = tvec(idx(1));
        clusters(ic).t_offset = tvec(idx(end));
        clusters(ic).mass = sum(tstat(idx));
        clusters(ic).mean_t = mean(tstat(idx));
    end
end

function T = build_gcp_table_electrode(S, events_all, fidx, tvec, win)
    rows_subj = [];
    rows_elec = [];
    rows_blk  = [];
    rows_ev   = {};
    rows_gcp  = [];
    tmask = tvec >= win(1) & tvec <= win(2);
    for ib = 1:numel(S)
        for ie = 1:numel(events_all)
            ev = events_all{ie};
            M = S(ib).meanTFREv.(ev);          % freq × time × elec
            g = squeeze(nanmean(nanmean(M(fidx, tmask, :), 1), 2)); % elec × 1
            g = g(:);
            nE = numel(g);
            rows_gcp  = [rows_gcp; g]; %#ok<AGROW>
            rows_subj = [rows_subj; S(ib).X_subj.(ev)(:)]; %#ok<AGROW>
            rows_elec = [rows_elec; S(ib).X_elec.(ev)(:)]; %#ok<AGROW>
            rows_blk  = [rows_blk; repmat(S(ib).blk, nE, 1)]; %#ok<AGROW>
            rows_ev   = [rows_ev; repmat({ev}, nE, 1)]; %#ok<AGROW>
        end
    end
    T = table(rows_subj, rows_elec, rows_blk, rows_ev, rows_gcp, ...
        'VariableNames', {'subj','elec','blk','event','GCP'});
end

function T = build_gcp_table_trial(S, events_all, fidx, tvec, win)
% Prefer explicit trial meta (X_subj_ch) if MAZE10 was re-run; else expand nEvents.
    rows_subj = [];
    rows_elec = [];
    rows_blk  = [];
    rows_ev   = {};
    rows_gcp  = [];
    tmask = tvec >= win(1) & tvec <= win(2);
    for ib = 1:numel(S)
        assert(~isempty(S(ib).ALL_TFR_MAT), 'ALL_TFR_MAT missing');
        for ie = 1:numel(events_all)
            ev = events_all{ie};
            A = S(ib).ALL_TFR_MAT.(ev);   % trial × freq × time
            if isfield(S(ib), 'X_subj_ch') && ~isempty(S(ib).X_subj_ch) && isfield(S(ib).X_subj_ch, ev)
                subj_tr = S(ib).X_subj_ch.(ev)(:);
                elec_tr = S(ib).X_elec_ch.(ev)(:);
                assert(numel(subj_tr) == size(A,1), 'X_subj_ch length mismatch for %s', ev);
            else
                nE = S(ib).nEvents.(ev)(:);
                subj_e = S(ib).X_subj.(ev)(:);
                elec_e = S(ib).X_elec.(ev)(:);
                assert(sum(nE) == size(A,1), ...
                    'nEvents sum (%d) != ALL_TFR_MAT rows (%d) for %s blk%d', ...
                    sum(nE), size(A,1), ev, S(ib).blk);
                subj_tr = [];
                elec_tr = [];
                for k = 1:numel(nE)
                    subj_tr = [subj_tr; repmat(subj_e(k), nE(k), 1)]; %#ok<AGROW>
                    elec_tr = [elec_tr; repmat(elec_e(k), nE(k), 1)]; %#ok<AGROW>
                end
            end
            g = squeeze(nanmean(nanmean(A(:, fidx, tmask), 2), 3)); % trial × 1
            g = g(:);
            rows_gcp  = [rows_gcp; g]; %#ok<AGROW>
            rows_subj = [rows_subj; subj_tr]; %#ok<AGROW>
            rows_elec = [rows_elec; elec_tr]; %#ok<AGROW>
            rows_blk  = [rows_blk; repmat(S(ib).blk, numel(g), 1)]; %#ok<AGROW>
            rows_ev   = [rows_ev; repmat({ev}, numel(g), 1)]; %#ok<AGROW>
        end
    end
    T = table(rows_subj, rows_elec, rows_blk, rows_ev, rows_gcp, ...
        'VariableNames', {'subj','elec','blk','event','GCP'});
end

function p_h = holm_correct(p)
    p = p(:)';
    m = numel(p);
    p_h = nan(size(p));
    [ps, ord] = sort(p);
    adj = nan(size(ps));
    for i = 1:m
        adj(i) = min(1, ps(i) * (m - i + 1));
    end
    % enforce monotonicity
    for i = m-1:-1:1
        adj(i) = min(adj(i), adj(i+1));
    end
    p_h(ord) = adj;
end

function plot_event_violins(T2, ev_order)
    hold on;
    cols = lines(numel(ev_order));
    for ie = 1:numel(ev_order)
        x = T2.GCP(T2.event == ev_order{ie} | strcmp(cellstr(T2.event), ev_order{ie}));
        % categorical safety
        if isempty(x)
            x = T2.GCP(strcmp(string(T2.event), ev_order{ie}));
        end
        swarmchart(ie*ones(size(x)), x, 12, cols(ie,:), 'filled', 'MarkerFaceAlpha', 0.35);
        plot(ie, nanmean(x), 'kd', 'MarkerFaceColor','k', 'MarkerSize', 8);
    end
    xlim([0.5 numel(ev_order)+0.5]);
    set(gca, 'XTick', 1:numel(ev_order), 'XTickLabel', ev_order);
    ylabel('OFC GCP');
    plot(xlim, [0 0], 'k:');
end

%% This code examines OFC oscill power at the time SWR occurs in HPC.
%% by contrasting the SWR-locked LPF with non-SWR-locked LPF with matching time window size
% statistical approach:
% ** channel-level **
% take OFC LFP +/= 200ms from HPC SWR peaks
% generate 401ms-long non-swr segments from task duration LFP
% random selection of non-swr segs, same # as swr segs
% generate channel-level t-maps for true and permuted data
% ** group-level **
% run mixed-model regressions that
%   1) tests the t-value intercept
%   2) includes (chan|subj)
% on true t-maps and permuted t-maps
% Use regression results from permuted t-maps to determine
% statistically significant intercept and cluster size

clc; clear all;
%% define paths

addpath(genpath('C:\Users\mgeva\Documents\GitHub\external\eeglab2023.1'))
addpath(genpath('C:\Users\mgeva\Documents\GitHub\Kamin_iEEG_MAZE_codes'));
addpath(genpath('C:\Users\mgeva\Documents\GitHub\external\KK_useful\useful'));
addpath('C:\Users\mgeva\Documents\GitHub\external\fieldtrip-20230613\fieldtrip-20230613\external\brewermap')
maze_set_path;
savedir = [dirs.TF 'pow_hpc_ripple_ofc_n7_grp/'];
%%
subjects = {'s01', 's06', 's08', 's09', 's12', 's14', 's15', 's16', 's17', 's18'};

% Discarding channels with high noise level
badRipples = {'s01_IndvRipples_RHH1','s08_IndvRipples_EEGRHB1REF','s09_IndvRipples_EEGRHH3REF'...
    's09_IndvRipples_EEGRHH3REF','s09_IndvRipples_EEGRHB1REF',...
    's12_IndvRipples_EEGLHH1REF',...
    's14_IndvRipples_LHD1_0005'};

goodRipplesubjects = ...
    {'s01' ,'s06' ,'s08' , 's15', 's16'  , 's17', 's18'};

peri_ripple_twin = 0.2; % in sec
nperm = 500;
osc_source = 'ofc'; %'ofc' or 'GM'
pac_phase_freq = 'alpha';

%% group vars: for TF-cell analysis
grpt_periripple_pow = [];
epairvec = [];
subjvec = [];
perm_tcell = {};
%% group vars: for event-wise analysis
c_count = 1;
agg_data = {};
fb_theta = [4 8];
fb_alpha = [8 10];
fb_lgamma = [30 50];
fb_mgamma = [30 100] ;

%%
for iSub = 1: length(subjects)
    subj = subjects{iSub};
    display(['Looking at subj ' subj]);
    
    subjdir_prep =  [dirs.prepfile subj '/'];   if ~exist(subjdir_prep); mkdir(subjdir_prep); end
    fname_chans = [subj '_ch_info.mat']; % bad channel idx
    fname_eeg = [subj '_maze_wtrig_pp.mat']; % eeg
    
    
    mm = matfile(fullfile(subjdir_prep, fname_chans));
    bchans = mm.bchans;
    mm = matfile(fullfile(subjdir_prep, fname_eeg));
    tmp = mm.data;
    ch_labels = tmp.elec.label;
    clear tmp
    
    load([dirs.banal, 'sq_table/', subj '_sq_table.mat']);
    task_onset = sq_table.time(1);
    task_offset = sq_table.time(end)+1;
    
    % load ripple (ied & gn data)
    load([dirs.ripple, '/', subj, '/' subj '_ripples.mat']);
    
    % get hpc-mpfc chan pairs
    chans_hpc = get_roi_chans(subj, 'HPC');
    chans_hpc = chans_hpc(ismember(chans_hpc, fieldnames(ripples)));
    % Remove HPC channels with faulty ripples (dominated by IED)
    rmv_ch= zeros(1,length(chans_hpc));
    for ii = 1:length(chans_hpc)
        for jj = 1:length(badRipples)
            
            if sum( (strfind(badRipples{jj}, chans_hpc{ii})) & ...
                    (strfind(badRipples{jj}, subj) ) )
                rmv_ch(ii) = 1;
            end
            
        end
    end
    chans_hpc = chans_hpc(~rmv_ch);
    
    chans_ofc = get_roi_chans(subj, 'OFC');
    coi_cidx = cell2mat(cellfun(@(x) find(ismember(ch_labels, x)), chans_ofc, 'UniformOutput', false));
    coi_cidx = setdiff(coi_cidx, bchans);
    if length(coi_cidx) < length(chans_ofc)
        error('bad channel involved')
    end
    
    % chans_ofc = chans_ofc(ismember(chans_ofc, fieldnames(ripples))); MGS
    % - why would we look for ripples in OFC?
    
    if ~isempty(chans_hpc) && ~isempty(chans_ofc) % only if there's at least one inter-regional elec pair
        for iHPC = 1: length(chans_hpc)
            % ripple-peaks
            ripple_peaks = ripples.(chans_hpc{iHPC}).peak;
            ripple_peaks = ripple_peaks(ripple_peaks>(task_onset-peri_ripple_twin) & ripple_peaks <(task_offset+peri_ripple_twin));
            
            %% get dpnts to reject
            rej = exc.(chans_hpc{iHPC});
            rej_HPC = [rej.GN rej.IED];
            
            
            if strcmp(osc_source, 'GM'), chans_ofc = 1; end
            for iOFC = 1: length(chans_ofc)
                clear rej
                
                % load OFC LFP data
                switch osc_source
                    case 'ofc'
                        load(fullfile(dirs.fourier, [subj '_' (chans_ofc{iOFC}) '.mat']));
                        display(['oscillatory activity in ' chans_ofc{iOFC} ' time-locked to ripples in ' chans_hpc{iHPC}]);
                        exc_ofc = load([dirs.ripple, '/', subj, '/' subj '_OF_IED.mat'],'exc');  % loading to get rej
                        exc_ofc = exc_ofc.exc;
                        
                        rej_OFC = exc_ofc.(chans_ofc{iOFC}); % Merge HPC and OFC rejected epochs
                        rej = [rej_HPC, rej_OFC.GN, rej_OFC.IED];
                        
                    case 'GM'
                        load(fullfile(dirs.fourier, [subj  '_GM.mat']));
                        rej = rej_HPC;
                        
                        fileList = dir(fullfile(dirs.fourier, [subj  '_*.mat']));
                        mm = matfile(fullfile(dirs.fourier, fileList(2).name));
                        hp = mm.hp;
                end
                
                %% get dpnts to reject [need hp.srate therefore under iOFC loop]
                % buffer noise width
                nbuffwin = 0.15;
                nbuffwin = round(hp.srate*nbuffwin);
                rej_buff = nan(length(rej), nbuffwin*2);
                for iR = 1: length(rej)
                    rej_buff(iR, :) = [rej(iR)-nbuffwin+1: rej(iR)+nbuffwin];
                end
                rej=unique(rej_buff);
                
                %% get twins around ripple peaks
                ripple_twin = [ripple_peaks-peri_ripple_twin ripple_peaks+peri_ripple_twin];
                ripple_swin = cat(2, dsearchn(lfp_time', ripple_twin(:, 1)), dsearchn(lfp_time', ripple_twin(:, 2)));
                nsamp = ripple_swin(1,2) - ripple_swin(1,1)+1;
                
                %% get event info that corresponds to each ripple peaks: for event-wise analysis
                var_etype = cell(length(ripple_peaks), 1);
                var_trl = nan(length(ripple_peaks), 1);
                var_block = nan(length(ripple_peaks), 1); % will repmat below
                var_impv = nan(length(ripple_peaks), 1);
                for ii = 1: length(ripple_peaks)
                    if sum(sq_table.time < ripple_peaks(ii)) ~= 0
                        sq_idx = find(sq_table.time == max(sq_table.time(sq_table.time < ripple_peaks(ii))));
                        var_etype{ii, 1} = sq_table.square{sq_idx};
                        var_trl(ii, 1) = sq_table.tidx(sq_idx);
                        var_block(ii, 1) = sq_table.block(sq_idx);
                        var_impv(ii, 1) = sq_table.opt_path_impv12(sq_idx);
                    end
                end
                emt = find(isnan(var_trl));
                var_etype(emt) = [];
                var_block(emt) = [];
                var_impv(emt) = [];
                var_trl(emt) = [];
                ripple_swin(emt,:) = [];
                
                %% epoch around hpc ripple peaks
                periripple_pow = nan(size(fourier.pow, 1), nsamp, size(ripple_swin, 1));
                periripple_phase = nan(size(fourier.phase, 1), nsamp, size(ripple_swin, 1));
                for ii = 1: size(ripple_swin, 1)
                    periripple_pow (:, :, ii) = fourier.pow(:, ripple_swin(ii, 1): ripple_swin(ii, 2));
                    periripple_phase (:, :, ii) = fourier.phase(:, ripple_swin(ii, 1): ripple_swin(ii, 2));
                    if sum(arrayfun(@(x) ( x>ripple_swin(ii, 1) & x< ripple_swin(ii, 2) ),rej)) ~= 0
                        periripple_pow(:, ii) = nan(size(periripple_pow(:, ii)));
                        periripple_phase(:, ii) = nan(size(periripple_phase(:, ii)));
                        display(['excluding event number ' num2str(ii)]);
                    end
                end
                outlier =[];
                tmp = squeeze(nanmean(periripple_pow, 2));
                outlier = cat(2, outlier, ...
                    find(sum(tmp > nanmean(tmp, 2) + nanstd(tmp, 0, 2)*10, 1) > 0), ...
                    find(sum(tmp < nanmean(tmp, 2) - nanstd(tmp, 0, 2)*10, 1) > 0));
                periripple_pow(:, :, outlier) = [];
                periripple_phase(:, :, outlier) = [];
                var_etype(outlier) = [];
                var_block(outlier) =[];
                var_impv(outlier) =  [];
                var_trl(outlier) = [];
                
                %% get equal number of non-ripple segments (same length as
                % periripple segs), free of noise
                for ii = 1: size(ripple_swin, 1)-1
                    tmp = [ripple_swin(ii, 2)+ round(0.05*hp.srate) : ripple_swin(ii+1,1) - round(0.05*hp.srate)];
                    nseg = floor(length(tmp)/nsamp);
                    for iSeg = 1: nseg
                        tmpseg = tmp((iSeg-1)*nsamp+1: iSeg*nsamp);
                        if isempty(intersect(rej, tmpseg))
                            IRI{ii,iSeg} = tmpseg;
                        end
                    end
                end
                segidx = find(cellfun(@(x) ~isempty(x), IRI));
                IRI = {IRI{segidx}};
                IRI_trl = nan(length(IRI), 1);
                IRI_block = nan(length(IRI), 1);
                IRI_impv = nan(length(IRI), 1);
                IRI_etype = cell(length(IRI), 1);
                for ii = 1: length(IRI)
                    ctr = median(IRI{ii})/hp.srate; % non-ripple event center (~ ripple peak)
                    sq_idx = find(sq_table.time == max(sq_table.time(sq_table.time < ctr)));
                    IRI_trl(ii, 1) = sq_table.tidx(sq_idx);
                    IRI_block(ii, 1) = sq_table.block(sq_idx);
                    IRI_impv(ii, 1) = sq_table.opt_path_impv12(sq_idx);
                    IRI_etype{ii, 1} = sq_table.square{sq_idx};
                end
                
                go = 1;
                while go
                    etype = unique(var_etype);
                    all_seg_idx = [];
                    for ii = 1: length(etype)
                        tmp_len1 = sum(var_block == 1 & strcmp(var_etype, etype{ii}));
                        tmp_len2 = sum(var_block == 2 & strcmp(var_etype, etype{ii}));
                        if ~isempty(var_block ==3), tmp_len3 = sum(var_block == 3 & strcmp(var_etype, etype{ii})); end
                        
                        segidx_blk1 = find(IRI_block == 1 & strcmp(IRI_etype, etype{ii}));
                        segidx_blk2 = find(IRI_block == 2 & strcmp(IRI_etype, etype{ii}));
                        if ~isempty(var_block ==3), segidx_blk3 = find(IRI_block == 3 & strcmp(IRI_etype, etype{ii})); end
                        
                        % block1
                        if length(segidx_blk1) > tmp_len1
                            segidx_blk1 = randsample(segidx_blk1, tmp_len1);
                        elseif length(segidx_blk1) < tmp_len1
                            pop_ridx = find(var_block == 1 & strcmp(var_etype, etype{ii}));
                            keep_ridx = randsample(pop_ridx, length(segidx_blk1));
                            pop_ridx = setdiff(pop_ridx, keep_ridx);
                            
                            periripple_pow (:, :, pop_ridx) = [];
                            periripple_phase (:, :, pop_ridx) = [];
                            var_etype(pop_ridx) = [];
                            var_trl(pop_ridx) = [];
                            var_block(pop_ridx) = [];
                            var_impv(pop_ridx) = [];
                        end
                        
                        % block2
                        if length(segidx_blk2) > tmp_len2
                            segidx_blk2 = randsample(segidx_blk2, tmp_len2);
                        elseif length(segidx_blk2) < tmp_len2
                            pop_ridx = find(var_block == 2 & strcmp(var_etype, etype{ii}));
                            keep_ridx = randsample(pop_ridx, length(segidx_blk2));
                            pop_ridx = setdiff(pop_ridx, keep_ridx);
                            
                            periripple_pow (:, :, pop_ridx) = [];
                            periripple_phase (:, :, pop_ridx) = [];
                            var_etype(pop_ridx) = [];
                            var_trl(pop_ridx) = [];
                            var_block(pop_ridx) = [];
                            var_impv(pop_ridx) = [];
                        end
                        
                        % block3
                        if ~isempty(var_block ==3)
                            if length(segidx_blk3) > tmp_len3
                                segidx_blk3 = randsample(segidx_blk3, tmp_len3);
                            elseif length(segidx_blk3) < tmp_len3
                                pop_ridx = find(var_block == 3 & strcmp(var_etype, etype{ii}));
                                keep_ridx = randsample(pop_ridx, length(segidx_blk3));
                                pop_ridx = setdiff(pop_ridx, keep_ridx);
                                
                                periripple_pow (:, :, pop_ridx) = [];
                                periripple_phase (:, :, pop_ridx) = [];
                                var_etype(pop_ridx) = [];
                                var_trl(pop_ridx) = [];
                                var_block(pop_ridx) = [];
                                var_impv(pop_ridx) = [];
                            end
                        end
                        all_seg_idx = cat(1, all_seg_idx, segidx_blk1, segidx_blk2, segidx_blk3);
                    end
                    
                    nonripple_pow = nan(size(periripple_pow, 1), size(periripple_pow, 2), length(all_seg_idx));
                    nonripple_phase = nan(size(periripple_pow, 1), size(periripple_pow, 2), length(all_seg_idx));
                    for iSeg = 1: length(all_seg_idx)
                        nonripple_pow(:, :, iSeg) = fourier.pow(:, IRI{all_seg_idx(iSeg)});
                        nonripple_phase(:, :, iSeg) = fourier.phase(:, IRI{all_seg_idx(iSeg)});
                    end
                    
                    outlier =[];
                    tmp = squeeze(nanmean(nonripple_pow, 2));
                    outlier = cat(2, outlier, ...
                        find(sum(tmp > nanmean(tmp, 2) + nanstd(tmp, 0, 2)*10, 1) > 0), ...
                        find(sum(tmp < nanmean(tmp, 2) - nanstd(tmp, 0, 2)*10, 1) > 0));
                    nonripple_pow(:, :, outlier) = [];
                    
                    if isequal(size(nonripple_pow), size(periripple_pow)), go = 0;  end
                end % end while
                clear IRI
                %% update event-wise analysis variables with non-ripple events
                var_block = cat(1, var_block, IRI_block(all_seg_idx));
                var_trl = cat(1, var_trl, IRI_trl(all_seg_idx));
                var_impv = cat(1, var_impv, IRI_impv(all_seg_idx));
                var_etype = cat(1, var_etype, IRI_etype{all_seg_idx});
                
                %% get epair-specific TF-map of t-values (peri- vs non-ripple events)
                diffmat = bsxfun(@minus, periripple_pow, nonripple_pow);
                tmat = nanmean(diffmat, 3) ./ (nanstd(diffmat, 0, 3)/sqrt(size(diffmat, 3)));
                
                grpt_periripple_pow = cat(3, grpt_periripple_pow, tmat);
                epairvec = [epairvec, iHPC*10+iOFC];
                subjvec = [subjvec, iSub];
                
                %% make permtmat: shuffle ripple occurrence
                perm_tmap_stacked = [];
                perm_all = cat(3, periripple_pow, nonripple_pow);
                for iPerm = 1: nperm
                    perm_all = perm_all(:, :, randperm(size(perm_all, 3)));
                    perm_rs = randsample(size(perm_all, 3), size(perm_all, 3)/2);
                    perm_a = perm_all(:, :, perm_rs);
                    perm_b = perm_all(:, :, setdiff([1:size(perm_all, 3)], perm_rs));
                    permdiffmat = bsxfun(@minus, perm_a, perm_b);
                    permtmat = nanmean(permdiffmat, 3) ./ (nanstd(permdiffmat, 0, 3)/sqrt(size(permdiffmat, 3)));
                    perm_tmap_stacked = cat(3, perm_tmap_stacked, permtmat);
                end
                perm_tcell = cat(1, perm_tcell, perm_tmap_stacked);
                
                %% compute event-wise measures
                
                % get fidx
                fidx.theta = dsearchn(hp.freqs, fb_theta');
                fidx.alpha = dsearchn(hp.freqs, fb_alpha');
                fidx.lgamma = dsearchn(hp.freqs, fb_lgamma');
                fidx.mgamma = dsearchn(hp.freqs, fb_mgamma');
                
                % pow: periripple
                y_theta = squeeze(nanmean(nanmean(periripple_pow(fidx.theta(1): fidx.theta(2), :, :), 2), 1));
                y_alpha = squeeze(nanmean(nanmean(periripple_pow(fidx.alpha(1): fidx.alpha(2), :, :), 2), 1));
                y_lgamma = squeeze(nanmean(nanmean(periripple_pow(fidx.lgamma(1): fidx.lgamma(2), :, :), 2), 1));
                y_mgamma = squeeze(nanmean(nanmean(periripple_pow(fidx.mgamma(1): fidx.mgamma(2), :, :), 2), 1));
                
                % pow: nonriripple
                y_theta2 = squeeze(nanmean(nanmean(nonripple_pow(fidx.theta(1): fidx.theta(2), :, :), 2), 1));
                y_alpha2 = squeeze(nanmean(nanmean(nonripple_pow(fidx.alpha(1): fidx.alpha(2), :, :), 2), 1));
                y_lgamma2 = squeeze(nanmean(nanmean(nonripple_pow(fidx.lgamma(1): fidx.lgamma(2), :, :), 2), 1));
                y_mgamma2 = squeeze(nanmean(nanmean(nonripple_pow(fidx.mgamma(1): fidx.mgamma(2), :, :), 2), 1));
                
                % pac
                if strcmp(pac_phase_freq,'alpha')
                    pac_phase = periripple_phase(fidx.alpha(1): fidx.alpha(2),:,:);
                    pac_phase2 = nonripple_phase(fidx.alpha(1): fidx.alpha(2),:,:);
                elseif strcmp(pac_phase_freq,'theta')
                    pac_phase = periripple_phase(fidx.theta(1): fidx.theta(2),:,:);
                    pac_phase2 = nonripple_phase(fidx.theta(1): fidx.theta(2),:,:);
                end
                pac_pow = squeeze(mean(periripple_pow(fidx.mgamma(1): fidx.mgamma(2),:,:), 1));
                pac_vals = nan(size(pac_phase, 1), size(pac_phase,3));
                pac_pow2 = squeeze(mean(nonripple_pow(fidx.mgamma(1): fidx.mgamma(2),:,:), 1));
                pac_vals2 = nan(size(pac_phase2, 1), size(pac_phase2,3));
                
                for iLF = 1: size(pac_phase, 1)
                    tmp = squeeze(pac_phase(iLF,:, :));
                    tmp2 = squeeze(pac_phase2(iLF,:, :));
                    pac_vals(iLF, :) = abs(mean(pac_pow.*exp(1i*tmp)));
                    pac_vals2(iLF, :) = abs(mean(pac_pow2.*exp(1i*tmp2)));
                end
                y_pac = nanmean(pac_vals, 1)';
                y_pac2 = nanmean(pac_vals2, 1)';
                
                % concat & clear
                y_theta = cat(1, y_theta, y_theta2); clear y_theta2;
                y_alpha = cat(1, y_alpha, y_alpha2);  clear y_alpha2;
                y_lgamma = cat(1, y_lgamma, y_lgamma2); clear y_lgamma2;
                y_mgamma = cat(1, y_mgamma, y_mgamma2);  clear y_mgamma2;
                y_pac = cat(1, y_pac, y_pac2);  clear y_pac2;
                
                % add ripple_occurrence, epair, and subj variable
                var_ripple = cat(1, repmat([1], [length(y_pac)/2, 1]), repmat([0], [length(y_pac)/2, 1]));
                var_epair = repmat([iSub*100+iHPC*10+iOFC], [length(y_pac), 1]);
                var_subj = repmat([iSub], [length(y_pac), 1]);
                
                %%
                agg_data{c_count, 1} = table(y_theta, y_alpha, y_lgamma, y_mgamma, y_pac, var_ripple, var_etype, var_block, var_trl, var_impv, var_epair, var_subj, ...
                    'VariableNames', {'theta', 'alpha', 'lg', 'mg', 'pac', 'ripple_occur', 'event', 'blk', 'trl', 'impv', 'epair', 'subj'});
                c_count = c_count+1;
                clear var_* y_* periripple* nonripple*
                
            end
        end
    end
end

%% generate meta_table for even-wise analysis
meta_table = agg_data{1};
for cc = 2:1: c_count-1
    meta_table = cat(1, meta_table, agg_data{cc});
end
disp('n-pairs')
obs_pairs = unique(meta_table.epair)

% save
header = fieldnames(meta_table)';
dcell = table2cell(meta_table);
dcell = cat(1, {header{1:size(dcell, 2)}}, dcell);

switch osc_source
    case 'ofc'
        cell2csv(fullfile(savedir, ['measure_newperm_meta_table_',pac_phase_freq,'.csv']), dcell);
        save(fullfile(savedir, ['5measure_newperm_meta_table_',pac_phase_freq,'.m']), 'meta_table', 'peri_ripple_twin','subjects')
    case 'GM'
        cell2csv(fullfile(savedir, ['measure_newperm_meta_table_GM_',pac_phase_freq,'.csv']), dcell);
        save(fullfile(savedir, ['5measure_newperm_meta_table_GM_',pac_phase_freq,'.m']), 'meta_table', 'peri_ripple_twin','subjects')
end

%% grp-level TF maps (Mixed-effects Model Regression & Permutation
mmt = nan(size(grpt_periripple_pow, 1)-2, size(grpt_periripple_pow, 2));
mmp = nan(size(grpt_periripple_pow, 1)-2, size(grpt_periripple_pow, 2));
for ii = 1: size(mmt, 1)
    for jj = 1: size(mmt, 2)
        Data = table(squeeze(grpt_periripple_pow(ii, jj, :)), epairvec', subjvec', 'VariableNames', {'tval', 'epair', 'subj'});
        model = fitlme(Data, 'tval ~ 1 + (epair|subj)');
        [beta, betanames, stats]  = fixedEffects(model);
        mmt(ii, jj) = stats.tStat(ismember(betanames.Name, '(Intercept)'));
        mmp(ii, jj) = stats.pValue(ismember(betanames.Name, '(Intercept)'));
    end
end

perm_mmt = nan(size(grpt_periripple_pow, 1)-2, size(grpt_periripple_pow, 2), nperm);
perm_mmp = nan(size(grpt_periripple_pow, 1)-2, size(grpt_periripple_pow, 2), nperm);
for iPerm = 1:nperm
    display(['running perm #' num2str(iPerm)]);
    a = cellfun(@(x) x(:, :, iPerm), perm_tcell, 'UniformOutput', false);
    for ii = 1: size(mmt, 1)
        parfor jj = 1: size(mmt, 2)
            b = cell2mat(cellfun(@(x) x(ii, jj), a,  'UniformOutput', false));
            Data = table(b , epairvec', subjvec', 'VariableNames', {'tval', 'epair', 'subj'});
            model = fitlme(Data, 'tval ~ 1 + (epair|subj)');
            [beta, betanames, stats]  = fixedEffects(model);
            perm_mmt(ii, jj, iPerm) = stats.tStat(ismember(betanames.Name, '(Intercept)'));
            perm_mmp(ii, jj, iPerm) = stats.pValue(ismember(betanames.Name, '(Intercept)'));
        end
    end
end
tthres_map = prctile(perm_mmt, 95, 3);
pmap_thres = (mmt > prctile(perm_mmt, 97.5, 3) | mmt < prctile(perm_mmt, 2.5, 3));

csize = nan(nperm, 1);
for iPerm = 1:  nperm
    perm_pmap = perm_mmp(:, :, iPerm) < .05;
    clustinfo = bwconncomp(perm_pmap);
    clust_info = cellfun(@numel,clustinfo.PixelIdxList);
    csize(iPerm, 1) = max(clust_info);
end
%
% clustinfo = bwconncomp(pmap_thres);
% clust_info = cellfun(@numel,clustinfo.PixelIdxList);
% whichclusters2remove = find(clust_info<prctile(csize, 95));
% % remove clusters
% for i=1:length(whichclusters2remove)
%     pmap_thres(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
% end

% alternative cluster thresholding
pmap_thres = (mmp < 0.005);
clustinfo = bwconncomp(pmap_thres);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);
whichclusters2remove = find(clust_info<35);
% remove clusters
for i=1:length(whichclusters2remove)
    pmap_thres(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
end

figure;
args = {[-0.2:1/hp.srate:0.2], hp.freqs(1:47), mmt}; %
map_handle = pcolor(args{:}); colormap(jet)
set(map_handle, 'EdgeColor', 'none');
set(gca, 'YScale', 'log');
set(gca,'YTick',[2 4 7 10 20 30 60 100 199], 'FontName', 'Arial', 'FontSize', 16); % for cnt
shading interp;
set(gca,'Box', 'off'); hold on;
[C, h] = contour([-0.2:1/hp.srate:0.2], hp.freqs(1:47), pmap_thres, [0 1], 'k', 'LineWidth', 2);
caxis([-5 5]); ylim([3 200]);colormap('Parula');
xlabel('Time (ms)',  'FontName', 'Arial', 'FontSize', 16);
ylabel('Frequency (Hz)', 'FontName', 'Arial', 'FontSize', 16);
cb = colorbar;
set(get(cb,'title'),'String', 't','FontName', 'Arial', 'FontSize', 12);
cb.FontSize = 12;
set(gcf,'name','ripple-grp-level TF maps')
savedir_pop = 'E:\MAZE_1.0\analysis\population_manuscript';
saveas(gcf, fullfile(savedir_pop,get(gcf,'name')),'tiff')

switch osc_source
    case 'ofc'
        saveas(gcf, [savedir, 'grp_stat_perm_matched_',pac_phase_freq], 'png');
    case 'GM'
        saveas(gcf, [savedir, 'grp_stat_perm_matched_gm_',pac_phase_freq], 'png');
end

% nperm = 500;
% perm_mmb = nan(size(grpt_periripple_pow, 1), size(grpt_periripple_pow, 2));
% for iPerm = 1: nperm
%     perm_data = shuffle(grpt_periripple_pow, 3);
%     tic
%     for ii = 1: size(mmt, 1)
%         for jj = 1: size(mmt, 2)
%             Data = table(squeeze(perm_data(ii, jj, :)), epairvec', subjvec', 'VariableNames', {'zpow', 'epair', 'subj'});
%             model = fitlme(Data, 'zpow ~ 1 + (epair|subj)');
%             [beta, betanames, stats]  = fixedEffects(model);
%             perm_mmb(ii, jj) = stats.tStat(ismember(betanames.Name, '(Intercept)'));
%     %         srp(ii, jj) = signrank(squeeze(periripple_pow_grp(ii, jj, :)));
%         end
%     end
%     toc
% end

%%

% Data = meta_table;
% obs_pairs = unique(Data.epair);
% etypes = unique(Data.event);
%
% for iType = 1: length(etypes)
%     tData = [];
%     for iBlk = 1: 2
%         for ii = 1: length(obs_pairs)
%             tmp = Data(Data.epair == obs_pairs(ii) & Data.blk == iBlk & strcmp(Data.event, etypes{iType}), :);
%             diffmat_pac = bsxfun(@minus, tmp(tmp.ripple_occur == 1, :).pac, tmp(tmp.ripple_occur == 0, :).pac);
%             tvals_pac(ii) = nanmean(diffmat_pac, 1) ./ (nanstd(diffmat_pac, 0, 1)/sqrt(size(diffmat_pac, 1)));
%             diffmat_alpha = bsxfun(@minus, tmp(tmp.ripple_occur == 1, :).alpha, tmp(tmp.ripple_occur == 0, :).alpha);
%             tvals_alpha(ii) = nanmean(diffmat_alpha, 1) ./ (nanstd(diffmat_alpha, 0, 1)/sqrt(size(diffmat_alpha, 1)));
%         end
%         tData = cat(1, tData, cat(2, obs_pairs, floor(obs_pairs/100), tvals_pac', tvals_alpha', repmat([iBlk], [length(obs_pairs) 1])));
%         clear tvals
%     end
%     tData_table{iType} = table(tData(:, 1), tData(:, 2), tData(:, 3), tData(:, 4), tData(:, 5),'VariableNames', {'epair', 'subj', 't_pac', 't_alpha', 'block'});
% end
%
% % print out stats
% for iType = 1: length(etypes)
%     display (etypes{iType})
%
%     tData = tData_table{iType};
%     tData(tData.t_pac==inf| tData.t_pac==-inf ,:) = [];
%     tData(isnan(tData.t_pac), :) = [];
%     tData.t_pacnorm = tData.t_pac ./ tData.t_alpha;
%     model = fitlme(tData, ['t_pacnorm ~ block + (epair|subj)'])
% end


%% pac after controlling for alpha
clear p pp
for ii_a = 1:2
    if ii_a == 1
        mm = matfile(fullfile(savedir, ['5measure_newperm_meta_table_GM_',pac_phase_freq,'.m']));
        meta_table = mm.meta_table;
    else
        mm = matfile(fullfile(savedir, ['5measure_newperm_meta_table_',pac_phase_freq,'.m']));
        meta_table = mm.meta_table;
    end
    
    Data = meta_table;
    % Data = Data(ismember(Data.event,{'X','D','G','N'}),:);
    obs_pairs = unique(Data.epair);
    
    tData = [];
    clear tvals_pac tvals_alpha tvals_mg
    for iBlk = 1
        for ii = 1: length(obs_pairs)
            tmp = Data(Data.epair == obs_pairs(ii) & Data.blk == iBlk & Data.impv == true, :);
            
            % MayaGS - using minus is problematic because PAC values + Alpha power are
            % positive and the distribution is biassed
            diffmat_pac = bsxfun(@minus, tmp(tmp.ripple_occur == 1, :).pac, mean(tmp(tmp.ripple_occur == 0, :).pac));
            % diffmat_pac = tmp(tmp.ripple_occur == 1, :).pac./tmp(tmp.ripple_occur == 0, :).pac;
            tvals_pac(ii) = nanmean(diffmat_pac, 1) ./ (nanstd(diffmat_pac, 0, 1)/sqrt(size(diffmat_pac, 1)));
            
            diffmat_alpha = bsxfun(@minus, tmp(tmp.ripple_occur == 1, :).alpha, mean(tmp(tmp.ripple_occur == 0, :).alpha));
            % diffmat_alpha = tmp(tmp.ripple_occur == 1, :).alpha./tmp(tmp.ripple_occur == 0, :).alpha;
            tvals_alpha(ii) = nanmean(diffmat_alpha, 1) ./ (nanstd(diffmat_alpha, 0, 1)/sqrt(size(diffmat_alpha, 1)));
            
            diffmat_mg = bsxfun(@minus, tmp(tmp.ripple_occur == 1, :).mg, mean(tmp(tmp.ripple_occur == 0, :).mg));
            % diffmat_mg = tmp(tmp.ripple_occur == 1, :).alpha./tmp(tmp.ripple_occur == 0, :).mg;
            tvals_mg(ii) = nanmean(diffmat_mg, 1) ./ (nanstd(diffmat_mg, 0, 1)/sqrt(size(diffmat_mg, 1)));
        end
        tData = cat(1, tData, cat(2, obs_pairs, floor(obs_pairs/100), tvals_pac', tvals_alpha', tvals_mg', repmat([iBlk], [length(obs_pairs) 1])));
        clear tvals
    end
    tData_table = table(tData(:, 1), tData(:, 2), tData(:, 3), tData(:, 4), tData(:, 5),tData(:, 6),...
        'VariableNames', {'epair', 'subj', 't_pac', 't_alpha','t_mg', 'block'});
    
    dist{1,ii_a} = tData_table.t_pac(~isoutlier(tData_table.t_pac,"mean") & ~isinf(tData_table.t_pac));
    
    dist{2,ii_a} = tData_table.t_alpha(~isoutlier(tData_table.t_alpha,"mean") & ~isoutlier(tData_table.t_pac,"mean")...
        & ~isinf(tData_table.t_alpha));

    dist{3,ii_a} = tData_table.t_mg(~isoutlier(tData_table.t_mg,"mean") & ~isoutlier(tData_table.t_pac,"mean") ...
         & ~isinf(tData_table.t_mg) );
    
    
    
    for iiS = 1:3
        p(iiS,ii_a) = signrank(dist{iiS,ii_a});
    end
    
end
clear pp
for iiS = 1:3
    [pp(iiS) h] = ranksum(dist{iiS,1},dist{iiS,2});
end

%% 
newA4figure('')
x_range = [-10 10];
y_range = [-0.22 0.3];

cmap = brewermap(6,'Dark2');
for ii = 1:3
    ii_a = 2 ; % OFC

    if ii == 1
        axes('position',[0.1 0.1 0.2 0.1])
        id = 2; %alpha

    elseif ii == 2
        axes('position',[0.1 0.26 0.2 0.1])
        id = 3; %gamma
    else
        axes('position',[0.1 0.42 0.2 0.1])
        ii_a = 2 ; % OFC
        id = 1; %pac

    end
    A = [dist{id,ii_a}'];
    raincloud_plot(A, 'color', cmap(id,:),'box_on',1,'alpha',0.5) ;
    set(gca,'xlim',x_range,'fontsize',8,'xtick',[-10,0,10],'ylim',y_range,'ytick',[0 0.3])
    hold all
    plot([0 0],get(gca,'ylim'),'k')
end
savename = 'distribution_ofc_alpha_PAC_gamma_rippleLocked_impv_blk1';
saveas(gcf, [savedir, savename], 'tiff');

%%
% newA4figure('')
% axes('position',[0.1 0.1 0.1 0.12])
% distributionPlot(dist{1,2},'histOpt',0,'showMM',2)
% savename = 'maze11_tval_pac_dist';
% saveas(gcf, [savedir, savename], 'png');

for ii_a = 1:2
    newA4figure('')
    if ii_a == 1
        A = [dist{1,2}' dist{1,1}']; % PAC
        B = [dist{2,2}' dist{2,1}']; % ALPHA
        ylabel_str = 'Ripple-locked Alpha-power (tstat)';
        savename = 'scatterhist_gm_ofc_alpha_PAC';
    else
        A = [dist{1,2}' dist{1,1}']; % PAC
        B = [dist{3,2}' dist{3,1}']; % GAMMA
        ylabel_str = 'Ripple-locked Gamma-power (tstat)';
        savename = 'scatterhist_gm_ofc_gamma_PAC';
        
    end
    groups = [ones(1,length(dist{1,2})),2*ones(1,length(dist{1,1}))];
    cmap = brewermap(4,'Paired');
    clr(1,:) = cmap(2,:); clr(2,:) = cmap(4,:);
    h = scatterhist(A,B,'Group',groups,'Nbins',10,'Marker','.','markersize',30,'color',clr);
    xlabel('Ripple-locked PAC (tstat)')
    ylabel(ylabel_str)
    legend('OFC','GM')
    axis([-5 5 -6 6]);
    hold on;
    boxplot(h(2),A,groups,'orientation','horizontal',...
        'label',{'OFC','GM'},'color',clr);
    boxplot(h(3),B,groups,'orientation','horizontal',...
        'label', {'OFC','GM'},'color',clr);
    set(h(2:3),'XTickLabel','');
    view(h(3),[270,90]);  % Rotate the Y plot
    axis(h(1),'auto');  % Sync axes
    hold on
    plot(h(1), [0 0],get(gca,'ylim'),'k')
    plot(h(1), get(gca,'xlim'),[0 0],'k')
    axes(h(3))
    hold on
    plot(h(3), [0 0],get(gca,'ylim'),'k')
    axes(h(2))
    hold on
    plot(h(2), [0 0],get(gca,'ylim'),'k')
    
    % saveas(gcf, [savedir, savename], 'tiff');
end

%% print out stats
model00 = fitlme(tData_table, ['t_pac ~ t_mg + block + (epair|subj)'],'DummyVarCoding' ,'effects');
stats = anova(model00)
model = fitlme(tData_table, ['t_pac ~  block + (epair|subj)'],'DummyVarCoding' ,'effects');
stats = anova(model)

figure
model0 = fitlme(tData_table, ['t_pac ~ t_alpha + (epair|subj)'],'DummyVarCoding' ,'effects');
stats = anova(model0)

cmap = brewermap(10,'Dark2');
obs_subj = unique(tData_table.subj);
for ii = 1:length(obs_subj)
    hold all
    id = find(tData_table.subj == obs_subj(ii));
    plot(tData_table.t_alpha(id), tData_table.t_pac(id),'.','color',cmap(ii,:),'markersize',25)
    % plot(tData_table.t_alpha(id), tData_table.t_pac(id),'.','color','k','markersize',25)
end


xlabel(sprintf(' Alpha - P_{ripple}-P_{control} (tStat)'))
ylabel(sprintf(' PAC_{ripple}-PAC_{control} (tStat)'))
title(sprintf('LME - t_{pac} ~ t_{alpha} + epair|subj, p = %2.2e',stats.pValue(2)))
axis([-5 5 -5 5])
hold all
plot(get(gca,'xlim'),[0 0],'k')
plot([0 0],get(gca,'ylim'),'k')
savename = 'maze11_tval_pac_alpha_scatter';
saveas(gcf, [savedir, savename], 'png');


figure
model0 = fitlme(tData_table, ['t_pac ~ t_mg + (epair|subj)'],'DummyVarCoding' ,'effects');
stats = anova(model0)

cmap = brewermap(10,'Dark2');
obs_subj = unique(tData_table.subj);
for ii = 1:length(obs_subj)
    hold all
    id = find(tData_table.subj == obs_subj(ii));
    plot(tData_table.t_mg(id), tData_table.t_pac(id),'.','color',cmap(ii,:),'markersize',25)
end
xlabel(sprintf(' Alpha - P_{ripple}-P_{control} (tStat)'))
ylabel(sprintf(' PAC_{ripple}-PAC_{control} (tStat)'))
title(sprintf('LME - t_{pac} ~ t_{gamma} + epair|subj, p = %2.2e',stats.pValue(2)))
axis([-4 4 -4 4])
hold all
plot(get(gca,'xlim'),[0 0],'k')
plot([0 0],get(gca,'ylim'),'k')
savename = 'maze11_tval_pac_mG_scatter';
saveas(gcf, [savedir, savename], 'png');

%%
Data = meta_table;
Data = Data(~isnan(Data.pac) & ~isnan(Data.alpha),:);
formula = 'pac ~ 1 + (epair|subj)';
lme1 = fitlme(Data,formula);
formula = 'pac ~  ripple_occur + (ripple_occur-1|epair)+(epair|subj)';
lme2 = fitlme(Data,formula);
results = compare(lme1,lme2);
disp(results.pValue)

%% copied from MAZE11_hpc_mpfc_pac_glm.m
% Comparing PAC modulation for improved and degraded trials

mm = matfile(fullfile(savedir, ['5measure_newperm_meta_table_',pac_phase_freq,'.m']));
meta_table = mm.meta_table;

data = meta_table;
data = data(~isnan(data.pac) & ~isnan(data.alpha),:);
model_ml = fitlme(data, ['impv ~ pac + event + blk + ripple_occur +  (epair|subj) '],'DummyVarCoding','effects');
anova(model_ml)


impv = 1;
data = meta_table;
data = data(data.impv == impv ,:);
model_ml0 = fitlme(data, ['pac ~ event + blk + (epair|subj) '],'DummyVarCoding','effects');
model_ml = fitlme(data, ['pac ~ event + blk  + ripple_occur +  (epair|subj) '],'DummyVarCoding','effects');
compare(model_ml0, model_ml)
anova(model_ml)

data = data(data.ripple_occur == 1 ,:);
data = data(strcmp(data.event,'N') ,:);
model_ml = fitlme(data, ['pac ~ blk  + (epair|subj) '],'DummyVarCoding','effects');
anova(model_ml)

data = meta_table;
impv = 0;
data = data(data.impv == impv ,:);
model_ml0 = fitlme(data, ['pac ~ event + blk + (epair|subj) '],'DummyVarCoding','effects');
model_ml = fitlme(data, ['pac ~ event + blk  + ripple_occur +  (epair|subj) '],'DummyVarCoding','effects');
compare(model_ml0, model_ml)
anova(model_ml)


data = meta_table;
data = data(~isnan(data.pac) & ~isnan(data.alpha),:);
data = data(data.ripple_occur == 1 ,:);
model_ml0 = fitlme(data, ['impv ~ blk + event +(epair|subj) '],'DummyVarCoding','effects');
model_ml = fitlme(data, ['impv ~ pac + blk + event +  (epair|subj) '],'DummyVarCoding','effects');
compare(model_ml0, model_ml)
anova(model_ml)

mm_gm = matfile(fullfile(savedir, ['5measure_newperm_meta_table_gm.m']))
meta_table_gm = mm_gm.meta_table;

%%
fig_name = sprintf('PAC_event_impv_degraded');

data = meta_table;
data = data(~isnan(data.pac) & ~isnan(data.alpha),:);
% data.pac = zscore(data.pac);
% impv = 1;
% data = data(data.impv == impv ,:);

newA4figure('')
EvTypePlot =  {'X','preG','N'};
x_w = 0.1;
y_w = 0.15;
positions(1,:) = [0.1 0.1 x_w y_w];
positions(2,:) = [0.22 0.1 x_w y_w];
positions(3,:) = [0.34 0.1 x_w y_w];
positions(4,:) = [0.48 0.1 x_w y_w];
positions(5,:) = [0.62 0.1 x_w y_w];
clear A A_se p
for ii = 1:length(EvTypePlot)
    axes('position',positions(ii,:))
    hold all
    title(EvTypePlot{ii})
    clear A B
    jj = 1;
    vec = data.pac(data.blk == 1 & data.impv == 1 & data.ripple_occur==1 & ~isnan(data.pac) & strcmp(data.event,EvTypePlot{ii}));
    vec1 = vec;
    A(jj) = nanmean(vec);
    A_se(jj) = nanstd(vec)/sqrt(length(vec));
    vec = data.pac(data.blk == 1  & data.impv == 0  & data.ripple_occur==1 & ~isnan(data.pac) & strcmp(data.event,EvTypePlot{ii}));
    vec2 = vec;
    
    jj = jj + 1;
    A(jj) = nanmean(vec);
    A_se(jj) = nanstd(vec)/sqrt(length(vec));
    %     vec = data.pac(data.ripple_occur==0 & ~isnan(data.pac) & strcmp(data.event,EvType{ii}) & data.blk == 2 ,:);
    %     jj = 1;
    %     B(jj) = nanmean(vec);
    %     B_se(jj) = nanstd(vec)/sqrt(length(vec));
    %     vec = data.pac(data.ripple_occur==1 & ~isnan(data.pac) & strcmp(data.event,EvType{ii}) & data.blk == 2 ,:);
    %     jj = jj + 1;
    %     B(jj) = nanmean(vec);
    %     B_se(jj) = nanstd(vec)/sqrt(length(vec));
    
    bb = bar(A);
    set(gca,'xtick',[1 2],'XTickLabel',{'Impv','Degraded'},'XTickLabelRotation',30 ,'ylim',[0 2],'YTick',[0 1])
    errorbar(1:2,A,A_se,'k.')
    % errorbar([2-0.15,2+.15],B,B_se,'k.')
    [p(ii),h] = ranksum(vec1,vec2);
    if p(ii) < 0.05
        plot(1.5,1.7,'k*','markersize',7)
    end
    if ii == 1;
        ylabel('PAC (Average)')
    end
end
% saveas(gcf, fullfile(savedir,fig_name),'tiff')
% close(gcf)

%%
fig_name = sprintf('PAC_event_rippleLocked_impv_blocks12');

data = meta_table;
data = data(~isnan(data.pac) & ~isnan(data.alpha),:);
% data.pac = zscore(data.pac);
impv = 1;
data = data(data.impv == impv ,:);

newA4figure('')
EvTypePlot =  {'X','preG','N'};
x_w = 0.1;
y_w = 0.15;
positions(1,:) = [0.1 0.1 x_w y_w];
positions(2,:) = [0.22 0.1 x_w y_w];
positions(3,:) = [0.34 0.1 x_w y_w];
positions(4,:) = [0.48 0.1 x_w y_w];
positions(5,:) = [0.62 0.1 x_w y_w];
clear A A_se p
for ii = 1:length(EvTypePlot)
    axes('position',positions(ii,:))
    hold all
    title(EvTypePlot{ii})
    clear A B
    jj = 1;
    vec = data.pac(data.blk ==1 & data.ripple_occur==1 & ~isnan(data.pac) & strcmp(data.event,EvTypePlot{ii}));
    vec1 = vec;
    A(jj) = nanmean(vec);
    A_se(jj) = nanstd(vec)/sqrt(length(vec));
    vec = data.pac(data.blk ==2   & data.ripple_occur==1 & ~isnan(data.pac) & strcmp(data.event,EvTypePlot{ii}));
    vec2 = vec;
    
    jj = jj + 1;
    A(jj) = nanmean(vec);
    A_se(jj) = nanstd(vec)/sqrt(length(vec));
    %     vec = data.pac(data.ripple_occur==0 & ~isnan(data.pac) & strcmp(data.event,EvType{ii}) & data.blk == 2 ,:);
    %     jj = 1;
    %     B(jj) = nanmean(vec);
    %     B_se(jj) = nanstd(vec)/sqrt(length(vec));
    %     vec = data.pac(data.ripple_occur==1 & ~isnan(data.pac) & strcmp(data.event,EvType{ii}) & data.blk == 2 ,:);
    %     jj = jj + 1;
    %     B(jj) = nanmean(vec);
    %     B_se(jj) = nanstd(vec)/sqrt(length(vec));
    
    bb = bar(A);
    set(gca,'xtick',[1 2],'XTickLabel',{'blk 1','blk 2'},'XTickLabelRotation',30 ,'ylim',[0 1.7])
    errorbar(1:2,A,A_se,'k.')
    % errorbar([2-0.15,2+.15],B,B_se,'k.')
    [p(ii),h] = ranksum(vec1,vec2);
    if p(ii) < 0.05
        plot(1.5,1.5,'k*','markersize',10)
    end
    if ii == 1;
        ylabel('PAC (Average)')
    end
end
% saveas(gcf, fullfile(savedir,fig_name),'tiff')
% close(gcf)

%% PAC is associated with locations in the maze Retrieval processes

Data_b2_impv = data(~isnan(data.pac) & ~isnan(data.alpha) & data.blk == 2 & data.impv == 1,:);
formula = 'pac ~ event + (epair|subj)';
lme1 = fitlme(Data_b2_impv,formula);
[p tbl stats] = anovan(Data_b2_impv.pac(:), {Data_b2_impv.event});
figure; multcompare(stats)

dataForDisplay = data(ismember(data.event,{'X','O','N','preG'}),:);
newA4figure('')
axes('position',[0.1 0.1 0.2 0.12])
distributionPlot(dataForDisplay.pac,'groups',dataForDisplay.event,'xValues',[1 3 2 4],'histOpt',0,'showMM',2)
axis([0 5 0 5])

fig_name = 'OFC_PAC_events';
saveas(gcf, fullfile(savedir,fig_name),'tiff')

%%

Data_ZS = data( ~isnan(data.pac) & ~isnan(data.alpha) ,:) ;
Data_ZS.pac = zscore(Data_ZS.pac);
Data_ZS.alpha = zscore(Data_ZS.alpha);
Data_ZS.mg = zscore(Data_ZS.mg);
Data_ZS.lg = zscore(Data_ZS.lg);
Data_ZS.theta = zscore(Data_ZS.theta);
formula = 'pac ~ event + (epair|subj)';
lme1 = fitlme(Data_ZS,formula);
formula = 'pac ~ ripple_occur + event  + (epair|subj)';
lme2 = fitlme(Data_ZS,formula);
results = compare(lme1,lme2);
disp(results.pValue)
[p tbl stats] = anovan(Data_ZS.pac(:), {Data_ZS.event});
figure; multcompare(stats)

Data_ZS = Data(~isnan(Data.pac) & ~isnan(Data.alpha) & ~isnan(Data.theta) & Data.impv == 1 ,:) ;
Data_ZS.pac = zscore(Data_ZS.pac);
Data_ZS.alpha = zscore(Data_ZS.alpha);
Data_ZS.mg = zscore(Data_ZS.mg);
Data_ZS.lg = zscore(Data_ZS.lg);
Data_ZS.theta = zscore(Data_ZS.theta);

formula = 'alpha ~ event * blk + (epair|subj)';
lme1 = fitlme(Data_ZS,formula);
formula = 'alpha ~ ripple_occur * event * blk  + (epair|subj)';
lme2 = fitlme(Data_ZS,formula);
results = compare(lme1,lme2);
disp(results.pValue)

[p tbl stats] = anovan(Data_ZS.lg(:), {Data_ZS.ripple_occur});
[p tbl stats] = anovan(Data_ZS.mg(:), {Data_ZS.ripple_occur});
[p tbl stats] = anovan(Data_ZS.alpha(:), {Data_ZS.event});
[p tbl stats] = anovan(Data_ZS.pac(:), {Data_ZS.event});

% Plot bar graph of these stats
maze11_plotStats()


xsz = 0.15;
ysz = 0.1;
position(1,:) = [0.1 0.1 xsz ysz];
position(2,:) = [0.1 position(1,2)+ysz*1.3 xsz ysz];
position(3,:) = [0.1 position(2,2)+ysz*1.3 xsz ysz];
position(4,:) = [0.1 position(3,2)+ysz*1.3 xsz ysz];
position(5,:) = [0.1 position(4,2)+ysz*1.3 xsz ysz];
position(6,:) = [0.1 position(5,2)+ysz*1.3 xsz ysz];

cmap = brewermap(6,'Dark2');
for ii = 1: 4
    axes('position',position(ii,:))
    raincloud_plot(A{ii}, 'color', cmap(ii,:),'box_on',1,'alpha',0.5) ;
    [p,h] = signrank(pac_sc);
    
    title( sprintf('ii: %d',ii))
    set(gca,'xlim',[-1 1])
    if ii == 1;    xlabel('zscore'); end
end


figure;
hold on
xsz = 0.15;
ysz = 0.1;
position(1,:) = [0.1 0.1 xsz ysz];
position(2,:) = [0.1 position(1,2)+ysz*1.3 xsz ysz];
position(3,:) = [0.1 position(2,2)+ysz*1.3 xsz ysz];
position(4,:) = [0.1 position(3,2)+ysz*1.3 xsz ysz];
position(5,:) = [0.1 position(4,2)+ysz*1.3 xsz ysz];
position(6,:) = [0.1 position(5,2)+ysz*1.3 xsz ysz];

Data_ZS = Data;
rmvId = isnan(Data_ZS.mg) | isnan(Data_ZS.lg);
Data_ZS(rmvId,:) = [];

cmap = brewermap(6,'Dark2');
for iEvType = 1: 6
    axes('position',position(iEvType,:))
    mg_sc = Data_ZS.mg(strcmp(Data_ZS.event,etypes{iEvType}) & Data_ZS.ripple_occur == 0) ;
    raincloud_plot(mg_sc, 'color', cmap(iEvType,:),'box_on',1,'alpha',0.5) ;
    [p,h] = signrank(mg_sc);
    
    title( sprintf('mg: %s, P = %2.2e',etypes{iEvType},p))
    set(gca,'xlim',[-1 6])
end


%% Fig. 4 panel A
str = {'all','impv','degraded'};

for ii_a = 1:3
    if ii_a == 1
        Data = meta_table;
    elseif ii_a == 2
        Data = meta_table;
        Data = Data(Data.impv == 1,:);
    elseif ii_a == 3
        Data = meta_table;
        Data = Data(Data.impv == 0,:);
    end
    
    obs_pairs = unique(Data.epair);
    etypes = unique(Data.event);
    
    tData_table = [];
    for iType = 1: length(etypes) -3
        tData = [];
        for iBlk = 1:2
            tvals_pac = []; tvals_alpha = [];
            for ii = 1: length(obs_pairs)
                didx = find(Data.epair == obs_pairs(ii) & Data.blk == iBlk & strcmp(Data.event, etypes{iType}));
                didx_impv = Data(didx,:).impv == 1;
                if length(didx) > 5 % only run if there is at least 5 observations
                    tmp = Data(didx, :);
                    diffmat_pac = bsxfun(@minus, tmp(tmp.ripple_occur == 1, :).pac, tmp(tmp.ripple_occur == 0, :).pac);
                    % diffmat_pac = tmp(tmp.ripple_occur == 1, :).pac - nanmean(tmp(tmp.ripple_occur == 0, :).pac);
                    tvals_pac = nanmean(diffmat_pac, 1) ./ (nanstd(diffmat_pac, 0, 1)/sqrt(size(diffmat_pac, 1)));
                    % diffmat_alpha = tmp(tmp.ripple_occur == 1, :).alpha - nanmean(tmp(tmp.ripple_occur == 0, :).alpha);
                    diffmat_alpha = bsxfun(@minus, tmp(tmp.ripple_occur == 1, :).alpha, tmp(tmp.ripple_occur == 0, :).alpha);
                    tvals_alpha = nanmean(diffmat_alpha, 1) ./ (nanstd(diffmat_alpha, 0, 1)/sqrt(size(diffmat_alpha, 1)));
                    tData = cat(1, tData, cat(2, obs_pairs(ii), floor(obs_pairs(ii)/100), tvals_pac, tvals_alpha, iType, iBlk));
                end
                clear tvals*
            end
        end
        tData_table = cat(1, tData_table, table(tData(:, 1), tData(:, 2), tData(:, 3), tData(:, 4), tData(:, 5),tData(:, 6),'VariableNames', {'epair', 'subj', 't_pac', 't_alpha', 'etype', 'block'}));
    end
    
    id = find(ismember(tData_table.etype, [3,4,6]));
    tData_table_mem = tData_table(id,:);
    model_mem = fitlme(tData_table_mem, 't_pac ~ block + etype + t_alpha + (epair|subj)');
    stats = anova(model_mem);
    pval(ii_a) = stats.pValue(3);
    
    % print out stats
    model0 = fitlme(tData_table, 't_pac ~ 1 + (epair|subj)')
    model1 = fitlme(tData_table, 't_pac ~ t_alpha + (epair|subj)')
    model = fitlme(tData_table, 't_pac ~ t_alpha + etype + block + (epair|subj)')
    compare(model0, model)
    
    cmap = brewermap(9,'Dark2');
    newA4figure('')
    pos(1,:) = [0.1 0.3 0.2 0.13];
    pos(2,:) = [0.1 0.1 0.2 0.13];
    clear ax
    for iblk = 1:2
        pac_m = nan(1, length(etypes)); alpha_m = nan(1, length(etypes));
        tData_table_blk = tData_table(tData_table.block == iblk,:);
        for iEvType = 1: length(etypes)
            pac_m(1, iEvType) = nanmean(tData_table_blk.t_pac (tData_table_blk.etype==iEvType, :));
            alpha_m(1, iEvType) = nanmean(tData_table_blk.t_alpha (tData_table_blk.etype==iEvType, :));
        end
        ax(iblk) = axes('position',pos(iblk,:));
        hold on
        for ii = [3:4,6]
            plot(alpha_m(ii),pac_m(ii),'.','color',cmap(ii,:),'markersize',25)
        end
        axis equal
        legend(etypes{[3:4,6]},'Location','EastOutside')
        
        xlabel('alpha diff(tStat)')
        ylabel('pac diff(tStat)')
        linkaxes(ax)
    end
    for iblk = 1:2
        axes(ax(iblk))
        plot([0 0], get(gca,'ylim'),'k','DisplayName','')
        plot( get(gca,'xlim'),[0 0], 'k','DisplayName','')
        title(sprintf('block %d',iblk))
    end
    suptitle(sprintf('changes between blocks with learning = %s',str{ii_a}));
    
    
    savedir = 'E:\MAZE_1.0\analysis\TF\pow_hpc_ripple_ofc_n7_grp\';
    savename = sprintf('tstat_pac_alpha_%s',str{ii_a});
    set(gcf,'name',savename)
    % saveas(gcf, [savedir, savename], 'tiff');
    
end




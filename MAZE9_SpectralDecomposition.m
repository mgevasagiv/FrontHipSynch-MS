%% run spectral decomposition
clear all; clc; 
%% all things PATH
restoredefaultpath; 
%% define paths
addpath(genpath('C:\Users\mgeva\Documents\GitHub\external\eeglab2023.1'))
addpath(genpath('C:\Users\mgeva\Documents\GitHub\Kamin_iEEG_MAZE_codes'));
addpath(genpath('C:\Users\mgeva\Documents\GitHub\external\KK_useful\useful'));
addpath('C:\Users\mgeva\Documents\GitHub\external\fieldtrip-20230613\fieldtrip-20230613\external\brewermap')
maze_set_path;
%%
subjects = {'s01', 's06', 's08', 's09', 's12', 's14', 's15', 's16', 's17', 's18'}; 
% subjects = {'s18'}; 
f_min = 2; f_max = 249; % 249 is the max (Nyquist)
rois = {'HPC'}; 
plot_photodiode = 1; 
hilbert_roichans = 0; 


for iSub = 1: length(subjects)
    subj = subjects{iSub}; 
    
    load(fullfile( [dirs.prepfile subj '/'], [subj '_ch_info.mat'])); % to get cmr_chans
    bchans = maze_get_bchans(subj); 
    bchans = cat(1, bchans.empty, bchans.pathological, bchans.ref);
    %% Channels
    if strcmp(rois{1}, 'GM')
        roichans = get_gm_chans (subj); 
        find_bchans = cell2mat(cellfun(@(x) sum(strcmp(x, bchans)), roichans, 'UniformOutput', false));         
        roichans = roichans(~find_bchans);
    else 
        for iRoi = 1: length(rois) 
            roichans{iRoi,:} = get_roi_chans(subj, rois{iRoi}); 
            exc = cell2mat(cellfun(@(x) ismember(x, bchans), roichans{iRoi,:}, 'UniformOutput', false));
            roichans{iRoi,:}(exc) = []; 
        end 
        if length(rois) == 2
            roichans =  cat(2, roichans{1}, roichans{2})
        elseif length(rois) == 1
            roichans =  roichans{1}
        end 
    end 
    %% Load pre-processed LFP data 
    dirs.prepfile = 'E:\MAZE_1.0\data\iEEG\preproc\';
    load([dirs.prepfile subj '\' subj '_maze_wtrig_pp.mat'])% cont.data, re-ref, downsampled, filtered 
%      load(['/Users/kaminkim/Documents/projects/iEEG_MAZE/data/iEEG/preproc/' subj '/' 'test.mat'])% cont.data, re-ref, downsampled, filtered 
    eeg_data = data.trial{1}; 
    %% Set frequency bands for spectral decomposiion 
    hp = get_hilbp (f_min, f_max, data); 
    %% Loop over roi chans 
    if hilbert_roichans
       if strcmp(rois{1}, 'GM')           
            gm_idx = setdiff(find(ismember(data.label, roichans)), cmr_chans);
            lfp = nanmean(eeg_data(gm_idx,:), 1); 
            [analytic_sig, fourier] = do_hilbert(lfp, hp);
            lfp_time = data.time{1}; 
            save([dirs.fourier subj '_GM'], 'analytic_sig', 'fourier', 'hp', 'lfp_time', '-v7.3');        
        else 
        for iChan = 1:length(roichans)
            cname = roichans{iChan}; 
            chanidx = strcmp(data.label, cname); 
            lfp = eeg_data(chanidx,:); 
            [analytic_sig, fourier] = do_hilbert(lfp, hp);
            lfp_time = data.time{1}; 
            save([dirs.fourier subj '_' cname], 'analytic_sig', 'fourier', 'hp', 'lfp_time', '-v7.3'); 
        end
       end 
    end 
    %% Do the same on photodiode and save plot 
    if plot_photodiode
        chanidx = find(strcmp(data.label, 'trigger')); 
        lfp = eeg_data(chanidx,:); 
%         EEG.data(trigger_idx, 214800:644000)
        [analytic_sig, fourier] = do_hilbert(lfp, hp);

        % epoch
        twin = [-0.2 0.2]; 
        swins = twin2samp (data, twin, subj); 
        for ii = 1: size(swins, 1)
            trial_pow(:, :, ii) = fourier.pow(:, swins(ii, 1): swins(ii, 2)); 
            trial_trace(ii, :) = lfp(:, swins(ii, 1): swins(ii, 2));    
        end 

        % plot 
        times = [twin(1)*1000: 1000/hp.srate: twin(2)*1000]; 
        figure; 
        subplot(1, 2, 1); plot(times, trial_trace'); 
        subplot(1, 2, 2);
    %     args = {times, hp.freqs(1:47,:), squeeze(mean(trial_pow(1:47, :, :), 3))}; %  
        args = {times, hp.freqs, squeeze(mean(trial_pow, 3))}; %  
        map_handle = pcolor(args{:}); colormap(jet)
        set(map_handle, 'EdgeColor', 'none'); 
        shading interp; colorbar
        saveas(gcf, fullfile(dirs.fourier, [subj, '_trigger' ]), 'png');  
    end 
    
    clear roichans data eeg_data hp analytic_sig fourier trial_*
end 

%% local functions
function [swin] = twin2samp (data, twin, subj)
    %% get trial onsets in sec
    [trigger_onsets] = get_trigger_onset (subj);

    win = twin*data.fsample; % in sample 
    sidx_onsets = dsearchn(data.time{1}', trigger_onsets); % data.time is in sec; finding sample indices at trial onsets
    swin = cat(2, sidx_onsets + win(1), sidx_onsets+win(2)); 
end 

function hp = get_hilbp (f_min, f_max, data)
    hp.bands = [];
    f = f_min;
    F = f_max; 
    while f < F
        hp.bands = cat(1,hp.bands,[f*0.95 f*1.05]);
        f = f*1.05/0.95; % for f-0.5f, f+0.5f
    end
    hp.freqs = mean(hp.bands,2);
    hp.srate = data.fsample;
    hp.order(1:size(hp.bands,1),1) = round((hp.srate./hp.bands(:,1))*3); % 3 times of the one cycle length for the lower bound frequency; MUST be in SAMPL
    % calc buffer length     
    hp.buffer_length_samp = round(max(hp.order)*1.5); % buffer length in sample 
end 

function [analytic_sig, fourier] = do_hilbert(lfp, hp)
    % initialize
    analytic_sig_buff = zeros(size(hp.bands,1), hp.buffer_length_samp*2+length(lfp));
    % buffer LFP trace with mirrored signal 
    lfp_buff = cat(2, ...
                cat(2, lfp(1) - flip(lfp(2: hp.buffer_length_samp+1), 2), lfp), ...
                lfp(end) - flip(lfp(end-hp.buffer_length_samp: end-1), 2)); %  1 x time(buffered) 
    % Hilbert transform
    for iFreq = 1:size(hp.bands,1)-1 % ~8-10s per freq
        [bcoeffs, acoeffs]= fir1(hp.order(iFreq), hp.bands(iFreq,:)/(hp.srate/2));
%         analytic_sig_buff(iFreq, :) = hilbert(filtfilt(bcoeffs, acoeffs, lfp_buff)')'; %  1 x time(buffered) 
        analytic_sig_buff(iFreq, :) = hilbert(filtfilt(bcoeffs, acoeffs, double(lfp_buff))')'; %  1 x time(buffered) 
    end 
    % truncate the buffer
    analytic_sig = analytic_sig_buff(:, hp.buffer_length_samp+1: end-hp.buffer_length_samp); 
    % extract power 
    fourier.pow = abs(analytic_sig).^2; 
    fourier.phase = angle(analytic_sig); 
end 
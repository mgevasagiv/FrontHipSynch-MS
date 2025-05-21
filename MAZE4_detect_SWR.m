clear all; clc;
%% all things PATH
restoredefaultpath;
addpath('C:\Users\mgeva\Documents\GitHub\external\eeglab2023.0')
addpath(genpath('C:\Users\mgeva\Documents\GitHub\Kamin_iEEG_MAZE_codes'));
maze_set_path;

%% select file -- get subject id and dir, fname
% select subject to proecess
subj_list = ...
    {'s01','s06','s08','s09','s12','s14','s15','s16','s17','s18'};

% cd(dirs.prepfile);
% dlist = dir('s*');
% [didx,~] = listdlg('PromptString', 'Select subject to process',...
%     'SelectionMode','single',...
%     'ListString', {dlist.name});

% get subject id and subjfolders
% subj = dlist(didx).name;

for ii_s = 1:length(subj_list)
    subj = subj_list{ii_s};
    subjdir_prep =  [dirs.prepfile subj '/'];   if ~exist(subjdir_prep); mkdir(subjdir_prep); end
    subjdir_ripple = [dirs.ripple subj '/'];    if ~exist(subjdir_ripple); mkdir(subjdir_ripple); end
    
    % file names
    fname_eeg = [subj '_maze_wtrig_pp.mat']; % eeg
    fname_chans = [subj '_ch_info.mat']; % bad channel idx
    fname_ripple=[subj '_ripples.mat'];
    
    % load files
    load(fullfile(subjdir_prep, fname_eeg));
    load(fullfile(subjdir_prep, fname_chans));
    
    % nan out bad channels
    data.trial{1}(bchans, :) = nan;
    [coi_cname, cref_cname] = maze_swr_coi(subj);
    coi_cidx = cell2mat(cellfun(@(x) find(ismember(data.elec.label, x)), coi_cname, 'UniformOutput', false));
    coi_cidx = setdiff(coi_cidx, bchans);
    cref_cidx = cell2mat(cellfun(@(x) find(ismember(data.elec.label, x)), cref_cname, 'UniformOutput', false));
    if length(cref_cidx) ~= length(coi_cidx)
        for kk = 1:length(cref_cidx)
            if (cref_cidx(kk)- coi_cidx) > 5
                rmv_idx(kk) = 1;
            else
                rmv_idx(kk) = 0;
            end
        end
        cref_cidx(logical(rmv_idx)) = [];
        clear rmv_idx
    end
    
    % run ripple detection
    [ripples,ripples_stat, exc] = detect_SWR(data, coi_cidx, cref_cidx); % this function preps data for SWR detection and calls the actual detection code (ripples_detection_excluding_IED_kk)
    save(fullfile(subjdir_ripple, fname_ripple), 'ripples', 'ripples_stat', 'exc');
    
    clear ripples ripple_stat exc
end
clear all; clc; 

ref_prep = 1; 
%% all things PATH
restoredefaultpath; 
addpath(genpath('/Users/kaminkim/Documents/codes/iEEG_MAZE/'));
addpath(genpath('/Users/kaminkim/Documents/codes/useful')); 
maze_set_path; 

%% select file -- get subject id and dir, fname
cd(dirs.prepfile); 
dlist = dir('s*'); 
[didx,~] = listdlg('PromptString', 'Select subject to process',...
                'SelectionMode','single',...
                'ListString', {dlist.name});
% get subject id            
subj = dlist(didx).name;
%% file names
subjdir =  [dirs.prepfile subj '/']; 
eegf_in = [subj '_maze_wtrig.set']; 
eegf_out = [subj '_maze_wtrig_pp.mat']; 
%% EEGLAB: prepare reref; load data
addpath('/Users/kaminkim/Documents/toolboxes/eeglab13_4_4b/');

% get cmr to be used for re-referencing 
if ref_prep
    prep_reref (subj, dirs, 'rf'); % separated out from this pipeline
end 
% load data and convert to ft format 
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename', eegf_in, 'filepath', subjdir);
%% FIELDTRIP: downsample, re-reference, and filter 
restoredefaultpath; % removing eeglab path before adding fieldtrip; rmpath doesn't work for removing eeglab from path :-/
addpath(genpath('/Users/kaminkim/Documents/toolboxes/fieldtrip/')); 
addpath(genpath('/Users/kaminkim/Documents/codes/iEEG_MAZE/'));
addpath(genpath('/Users/kaminkim/Documents/codes/useful')); 
maze_set_path; 

% convert data to ft format 
data = eeglab2fieldtrip(EEG, 'preprocessing', 'none'); clear EEG; 

% downsample 
cfg = [];
cfg.resamplefs = 512; 
data_ds = ft_resampledata(cfg, data); 

% rereference & filter 
load([dirs.prepfile, subj, '/', subj '_ch_info.mat']); 
cfg = [];
cfg.channel     =  {'all', '-trigger'}; 
cfg.dftfilter   = 'no'; % 
cfg.dftfreq     =[60 120 180 240]; 
cfg.bsfilter    = 'yes'; % US line noise bandstop frequency range. 
cfg.bsfreq     = [59 61; 119 121; 179 181; 239 241]; % US line noise bandstop frequency range. 
% cfg.lpfilter   = 'no'; % no need - downsampled to 500 Hz 
% cfg.hpfilter   = 'yes'; % high-pass in order to get rid of low-freq trends
% cfg.hpfiltord  = 3; %high pass filter order specified as 3
% cfg.hpfreq     = 1; %hp frequency 1. Pass anything <1 Hz 
% cfg.lpfiltwintype = 'hann'; %setting hamming window 
% cfg.hpfiltwintype = 'hann'; %setting up hamming window 
cfg.padding = 2; 

cfg.padtype = 'mirror'; 
cfg.reref         =  'yes'; 
cfg.refchannel    = data.elec.label(cmr_chans); 
cfg.refmethod     = 'avg'; 
data_ds_ref_filt     = ft_preprocessing(cfg,data_ds);
clear data; 
data = data_ds_ref_filt; 
data.trial{1} = cat(1, data.trial{1}, data_ds.trial{1}(find(cell2mat(cellfun(@(x) strcmp(x, 'trigger'), data_ds.label, 'UniformOutput', false))),:));
data.label{find(cell2mat(cellfun(@(x) strcmp(x, 'trigger'), data_ds.label, 'UniformOutput', false)))} = 'trigger'; 
% cfg = [];
% cfg.viewmode = 'vertical'; 
% ft_databrowser(cfg, data)
% save(fullfile(subjdir, 'pp_test_bson_nootherfilt_240_1000Hz'), 'data');

save(fullfile(subjdir, eegf_out), 'data',  '-v7.3');


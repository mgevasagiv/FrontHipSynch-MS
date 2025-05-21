%% Instead of the approach of epoching ripple rate time series,
% obtain event-related ripple rats using PSTH (centering bins around the
% the event of interest)
% Edited by MGS Dec 2023

clc; clear all;
%% define paths
addpath(genpath('C:\Users\mgeva\Documents\GitHub\external\eeglab2023.1'))
addpath(genpath('C:\Users\mgeva\Documents\GitHub\Kamin_iEEG_MAZE_codes'));
maze_set_path;
%% behavior-driven ripples (linking behav & eeg tstamps)
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
    % load EEG to get the trigger timestamps- timestamps are in ms; srate = raw
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = pop_loadset( 'filename', [subj '_maze_wtrig.set'], 'filepath', [dirs.prepfile subj '/']);
    EEG = eeg_checkset( EEG ); close all;
    
    % convert timestamps into s
    trigger_onsets = [EEG.event.latency]' ; % all NewMaze onsets, in SAMPLE
    trigger_onsets = trigger_onsets/EEG.srate; % now in s
    
    % load behavioral tagging, get onset time for each event type in s
    % note that 'tidx' in binfo is actually 'event idx'
    load([dirs.behavioral, '\binfo\', subj '_binfo.mat']); % in ms
    % load matching decition points for D and mD
    load([dirs.banal, '\nav_impv_v3\', subj '_matching_decpnts.mat']); %MGS - what version is the relevant one?
    
    events = fieldnames(square);
    for i = 1 : length(events)
        tmp = square.(events{i});
        if strcmp(events{i}, 'D') % remove non-matching decrrision pnts
            [tmp, ~] = func_reduce_D (tmp, dpmatch);
        end
        for iEvent = 1: length(tmp)
            thisEvent = (tmp(iEvent).block-1)*30 + tmp(iEvent).tidx;
            EventOnsets.(events{i})(iEvent) = trigger_onsets(thisEvent)+tmp(iEvent).tfromMazeOn/1000; % in s
        end
        clear tmp
    end
    
    tmpm = {'D', 'N', 'O', 'X', 'G'};
    for i = 1: length(tmpm)
        tmp = square.(tmpm{i});
        if strcmp(tmpm{i}, 'D') % remove non-matching decision pnts
            [tmp, ~] = func_reduce_D (tmp, dpmatch);
        end
        movestr = ['m' tmpm{i}];
        for iEvent = 1: length(tmp)
            thisEvent = (tmp(iEvent).block-1)*30 + tmp(iEvent).tidx;
            EventOnsets.(movestr)(iEvent) = trigger_onsets(thisEvent)+tmp(iEvent).tfromMazeOn/1000+tmp(iEvent).rt/1000; % in s
        end
    end
    
    % load ripple data
    load([dirs.ripple, '/', subj, '/' subj '_ripples.mat'])
    evlist = fieldnames(EventOnsets);
    for Eidx = 1: length(evlist)
        EOI = evlist{Eidx}; % event of interest
        fsize = 5; % triangular filter width (dpoints)
        twin = [1.5 1.5]; % in sec
        % how many bins?
        % bin_width shoule be about 100-150 samples, which is equivelant to a width of 200-300ms (Norman)
        % We may need a smaller window given the faster-tempoed nature of our task
        sws_s = 0.05;
        chans = fieldnames(ripples);
        for iElec = 1: length(chans)
            for iEvent = 1: length(EventOnsets.(EOI))
                edges(iEvent,:) = [EventOnsets.(EOI)(iEvent) - twin(1) : sws_s: EventOnsets.(EOI)(iEvent) + twin(2)];
                count(iEvent, :) = histcounts(ripples.(chans{iElec}).peak,  edges(iEvent,:));
                r = count(iEvent, :) / sws_s; % per sec
                EventRR.(chans{iElec})(iEvent, :) = nanfastsmooth(r, fsize-1, 2, 0.5); % 2 =triangular filter
                clear r;
            end
        end
        save ([dirs.ripplerate_ER, 'rr_' EOI '/', subj '_RippleRate_ER_' EOI '.mat'], 'EventRR', 'edges', 'count', 'twin', 'sws_s');
        clear EventRR;
    end
    
    evlist = {'X_1ext', 'X_2ext'};
    for Eidx = 1: length(evlist)
        EOI =  evlist{Eidx};
        tmp_onsets = EventOnsets.X([square.X.nexits] == Eidx);
        fsize = 5; % triangular filter width (dpoints)
        twin = [1.5 1.5]; % in sec
        sws_s = 0.05;
        chans = fieldnames(ripples);
        for iElec = 1: length(chans)
            for iEvent = 1: length(tmp_onsets)
                edges(iEvent,:) = [tmp_onsets(iEvent) - twin(1) : sws_s: tmp_onsets(iEvent) + twin(2)];
                count(iEvent, :) = histcounts(ripples.(chans{iElec}).peak,  edges(iEvent,:));
                r = count(iEvent, :) / sws_s; % per sec
                EventRR.(chans{iElec})(iEvent, :) = nanfastsmooth(r, fsize-1, 2, 0.5); % 2 =triangular filter
                clear r;
            end
        end
        save ([dirs.ripplerate_ER, 'rr_' EOI '/', subj '_RippleRate_ER_' EOI '.mat'], 'EventRR', 'edges', 'count', 'twin', 'sws_s');
        clear EventRR;
    end
    
end
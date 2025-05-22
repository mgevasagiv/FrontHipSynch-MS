% input: roi = 'HPC' or 'OFC'
%
% Generates TFR 
% output: sprintf('CenteredAvgEvents_allElectrodesSubj_%s_%d.mat',roi,cond)
%        sprintd('CenteredAvgEvents_allElectrodesSubj_%s_%d_TFRMAT.mat',roi,impvBlk2)
% 
% Conditions - 
% 0 - block 1 (improved trials)
% 1 - block 2 (improved trials)
% 3 - block 1 (degraded trials)
% 4 - block 2 (degraded trials)

function [] = MAZE10_calc_TFR_event(roi)

% maze_set_path;
github_path = 'C:\Users\mgeva\Documents\GitHub\';
addpath(genpath(fullfile(github_path, 'closedLoop-code')));
addpath(genpath(fullfile(github_path, 'external\fieldtrip-20230613\fieldtrip-20230613')))
rmpath(('C:\Users\mgeva\Documents\GitHub\external\fieldtrip-20230613\fieldtrip-20230613\compat\matlablt2010b'))
addpath('C:\Users\mgeva\Documents\GitHub\external\KK_useful\useful')
addpath('E:\Dropbox\Code\Useful_code')

maze_set_path;

%% Edit below this line to specify subjects included in analysis
subjects = {'s01', 's06', 's08', 's09', 's12', 's14', 's15', 's16', 's17', 's18'};
min_n_events = 5;

params.mgamma_range = [50,100]; % mgamma range

%% 

for iiCond = 3
    if iiCond == 1
        impvBlk2 = 1;
    elseif iiCond == 2
        impvBlk2 = 0;
    elseif iiCond == 3
        impvBlk2 = 3;
    elseif iiCond == 4
        impvBlk2 = 4;
    end
    
    savedir = [dirs.TF 'TFR_grp/' roi '/'];
    evlist = {'X', 'D', 'G', 'E', 'O', 'N', 'preG', 'preO','null'};
    
    
    %%
    srate= 512; % Hz - assuming this holds for all final datasets
    win = [0.5 0.5]; % in sec
    
    cnt = 0; % number of electrodes in the calculation
    
    nRip = [];
    clear avgRipple stdRipple
    clear meanTFRRip meanTFRRip spectRip nRip
    
    Y_full = [];
    X_ttype = [];
    X_blk = [];
    X_elec = [];
    X_subj = [];
    chCnt = zeros(1,length(evlist));
    
    for iEtype = 1:length(evlist)
        chCnt(iEtype) = 1;
        meanTFREv.(evlist{iEtype}) = [];
        X_blk.(evlist{iEtype}) = [];
        X_ttype.(evlist{iEtype}) = [];
        X_elec.(evlist{iEtype}) = [];
        X_subj.(evlist{iEtype}) = [];
        avgLFP.(evlist{iEtype})  = [];
        stdLFP.(evlist{iEtype}) = [];
        spectEv.(evlist{iEtype}) = [];
        nEvents.(evlist{iEtype}) = [];
        GammaBehavTable.(evlist{iEtype}) = [];
        
        GammaBand_trial_all.(evlist{iEtype}) = [];
        ALL_TFR_MAT.(evlist{iEtype}) = [];
        opt_path_impv12_trial.(evlist{iEtype}) = [];
        BehavD2D.(evlist{iEtype}) = [];
        BehavD2G.(evlist{iEtype}) = [];
        X_elec_ch.(evlist{iEtype})  = [];
        X_subj_ch.(evlist{iEtype}) = [];
        
    end
    
    for iSub = 1: length(subjects)
        
        subj = subjects{iSub};
        clist = get_roi_chans(subj, roi);
        if isempty(clist); continue; end
        
        % load behavioral and noise data
        load([dirs.banal, 'sq_table/', subj '_sq_table.mat']);
        if strcmp(roi, 'OFC')
            load([dirs.ripple, '/', subj, '/' subj '_OF_IED.mat']); clear ripples*; % loading to get rej
        end
        if strcmp(roi, 'HPC')
            load([dirs.ripple, '/', subj, '/' subj 'HPC_IED.mat']); clear ripples*; % loading to get rej
        end
        % load analytic signal channel names ('analytic_sig', 'fourier', 'hp', 'lfp_time'); lfp_time in sec
        clist = clist(ismember(clist, fieldnames(exc)')); % These fieldnames do not include bad channels
        
        task_onset = sq_table.time(1);
        task_offset = sq_table.time(end)+1;
        
        if impvBlk2 == 1
            sq_table = sq_table(sq_table.opt_path_impv12 == 1  & sq_table.block == 2 ,:);
        elseif impvBlk2 == 0
            sq_table = sq_table(sq_table.opt_path_impv12 == 1  & sq_table.block == 1 ,:);
        elseif impvBlk2 == 3
            sq_table = sq_table(sq_table.opt_path_impv12 == 0  & sq_table.block == 1 ,:);
        elseif impvBlk2 == 4
            sq_table = sq_table(sq_table.opt_path_impv12 == 0  & sq_table.block == 2 ,:);
        end
        
        dirs.prepfile = 'E:\MAZE_1.0\data\iEEG\preproc\';
        load([dirs.prepfile subj '\' subj '_maze_wtrig_pp.mat'])% cont.data, re-ref, downsampled, filtered
        eeg_data = data.trial{1};
        
        for iElec = 1: length(clist)
            % load analytic signal ('analytic_sig', 'fourier', 'hp', 'lfp_time'); lfp_time in sec
            % load(fullfile(dirs.fourier, [subj '_' clist{iElec} '.mat']));
            
            cname = clist{iElec};
            
            if iElec > 1
                if strcmp(cname(1:end-1),clist{iElec-1}(1:end-1))
                    continue
                end
            end
            chanidx = strcmp(data.label, cname);
            lfp = eeg_data(chanidx,:);
            
            % load noise data
            rej = exc.(clist{iElec});
            rej = [rej.GN rej.IED];
            % buffer noise width
            nbuffwin = 0.15;
            nbuffwin = round(srate*nbuffwin);
            rej_buff = nan(length(rej), nbuffwin*2);
            for iR = 1: length(rej)
                rej_buff(iR, :) = [rej(iR)-nbuffwin+1: rej(iR)+nbuffwin];
            end
            rej=unique(rej_buff);
            
            % set params
            t = double(data.time{1}); % in sec
            srate = data.fsample;
            if srate ~= 512
                error('data needs to be downsampled')
            end
            center = win(1)*data.fsample+1;
            
            params.freqRangeForAvgSpec = [1:250];
            params.timeBeforeAfterEventRipSpec = 0.75; % sec
            params.timeForBaselineRip = 1; % sec
            params.minNCycles = 5;
            params.freqoiForAvgSpec = [0:0.5:10];
            pacCalc = PACCalculator;
            pacCalc.samplingRate = srate;
            pacCalc.freqRange = params.freqRangeForAvgSpec;
            pacCalc.timeBeforeAfterEvent = params.timeBeforeAfterEventRipSpec; %seconds
            pacCalc.timeForBaseline = params.timeForBaselineRip; %seconds, from Starestina et al
            pacCalc.minNCycles = params.minNCycles;
            
            for iEtype = 1: length(evlist)
                
                
                %% get epoch window for this trial type
                if ~strcmp(evlist{iEtype}, 'E')
                    trl_twin = cat(2, sq_table.time(strcmp(sq_table.square, evlist{iEtype})), ...
                        sq_table.time(strcmp(sq_table.square, evlist{iEtype}))+ sq_table.rt(strcmp(sq_table.square, evlist{iEtype}))/1000);
                else
                    trl_twin = cat(2, sq_table.time(strcmp(sq_table.square, evlist{iEtype})), ...
                        sq_table.time(strcmp(sq_table.square, evlist{iEtype}))+ 1);
                end
                rowId = find(strcmp(sq_table.square, evlist{iEtype}));
                d2dec = sq_table.d2dec(rowId);
                d2goal = sq_table.d2goal(rowId);
                opt_path_impv12 = sq_table.opt_path_impv12(rowId);
                
                
                % trl_swin = cat(2, dsearchn(lfp_time', trl_twin(:, 1)), dsearchn(lfp_time', trl_twin(:, 2)));
                t1 = dsearchn(data.time{1}', trl_twin(:, 1));
                trl_swin = cat(2, t1, t1+srate*params.timeBeforeAfterEventRipSpec); % range after event
                buf_control = srate*(params.timeBeforeAfterEventRipSpec + params.timeForBaselineRip); % range before event, including control
                
                %% epoch
                peaks_local = []; rmv_id = [];
                for ii = 1: size(trl_swin, 1)
                    if sum(arrayfun(@(x) ( x>(trl_swin(ii, 1)-buf_control) & (x < trl_swin(ii, 2)) ),rej)) ~= 0
                        % SquaredurPow.(evlist{iEtype})(:, ii) = nan(size(SquaredurPow.(evlist{iEtype})(:, ii)));
                        rmv_id = [rmv_id, ii];
                        display(['excluding event number ' num2str(ii)]);
                        continue
                    else
                        peaks_local = [peaks_local trl_swin(ii,1)];
                    end
                end
                if length(peaks_local) < min_n_events continue; end
                
                [meanTFREv_ch, meanEpochs, stdEpochs, allEpochs, nEpochs] = pacCalc.plotAvgSpecDiff( lfp(:)', peaks_local);
                
                tvec = -params.timeBeforeAfterEventRipSpec:1/srate:params.timeBeforeAfterEventRipSpec; % sec
                switch evlist{iEtype}
                    case 'N'
                        time_to_test = [0.2 0.35];
                    case 'D'
                        time_to_test = [-0.35 -0.15];
                    case 'O'
                        time_to_test = [0.2 0.35];
                    case 'G'
                        time_to_test = [0.2 0.35];
                    case 'X'
                        time_to_test = [0.2 0.35];
                    case 'E'
                        time_to_test = [0.2 0.35];
                    case 'preG'
                        time_to_test = [0.2 0.35];
                    case 'preO'
                        time_to_test = [0.2 0.35];
                    case 'null'
                        time_to_test = [0.2 0.35];
                end
                id1_t_test = find(time_to_test(1)<tvec,1);
                id2_t_test = find(time_to_test(2)<tvec,1);
                id1_f_test = find(params.mgamma_range(1)<params.freqRangeForAvgSpec,1);
                id2_f_test = find(params.mgamma_range(2)<params.freqRangeForAvgSpec,1);
                
                clear GammaBand_trial
                for ii = 1:size(allEpochs,1)
                    GammaBand_trial(ii) = nanmean(nanmean( allEpochs(ii, id1_f_test:id2_f_test,id1_t_test:id2_t_test) ));
                end
                d2dec(rmv_id) = [];
                d2goal(rmv_id) = [];
                opt_path_impv12(rmv_id) = [];
                GammaBand_trial_all.(evlist{iEtype}) = cat(1, GammaBand_trial_all.(evlist{iEtype}), GammaBand_trial(:));
                opt_path_impv12_trial.(evlist{iEtype}) = cat(1, opt_path_impv12_trial.(evlist{iEtype}), opt_path_impv12(:));
                BehavD2D.(evlist{iEtype}) = cat(1, BehavD2D.(evlist{iEtype}), d2dec(:));
                BehavD2G.(evlist{iEtype}) = cat(1, BehavD2G.(evlist{iEtype}), d2goal(:));
                X_elec_ch.(evlist{iEtype}) = cat(2, X_elec_ch.(evlist{iEtype}), repmat(10*iSub+iElec,1,length(d2dec)));
                X_subj_ch.(evlist{iEtype}) = cat(2, X_subj_ch.(evlist{iEtype}), repmat(iSub,1,length(d2dec)));
                ALL_TFR_MAT.(evlist{iEtype}) = cat(1, ALL_TFR_MAT.(evlist{iEtype}), allEpochs);
                
                plotChSpectogram = 0;
                if (plotChSpectogram)
                    fname_fig = sprintf('eventLockedTFR_%s_%s_%s',subjects{iSub},cname,evlist{iEtype});
                    f = newA4figure(fname_fig)
                    set(gcf,'DefaultAxesFontSize',12);
                    
                    axes('position',[0.12,0.6,0.3,0.3])
                    % plot(avgWeighted, 'Color', 'g',  'LineWidth', 2); hold on;
                    tvec = -.75:1/srate:.75; % sec
                    shadedErrorBar(tvec,meanEpochs,stdEpochs/sqrt(nEpochs))
                    axis tight
                    %xlim([t(1) t(center*2-1)]);
                    xticks([tvec(1) tvec(end)]);
                    % axis off
                    hold on
                    plot([0 0], get(gca,'ylim'),'k')
                    %plot([t(1) t(1)+0.1],-10*ones(1,2),'linewidth',10,'color','k')
                    %plot([t(1) t(1)],[10 20],'linewidth',10,'color','k')
                    
                    axes('position',[0.12,0.15,0.35,0.3])
                    imagesc(tvec,pacCalc.freqRange,meanTFREv_ch)
                    hold on
                    plot([0 0], get(gca,'ylim'),'k')
                    
                    axis xy
                    colorbar
                    saveas(f,fullfile(savedir,fname_fig),'png')
                    close(f)
                end
                
                meanTFREv.(evlist{iEtype}) = cat(3, meanTFREv.(evlist{iEtype}) , meanTFREv_ch);
                
                % RawSig.(evlist{iEtype}) = cat(1, RawSig.(evlist{iEtype}) , LFP_chunks);
                
                % Y_full = cat(1, Y_full, SquaredurPow.(evlist{iEtype})');
                % X_blk.(evlist{iEtype}) = cat(1,  X_blk.(evlist{iEtype}),  sq_table.block(strcmp(sq_table.square, evlist{iEtype})));
                % X_ttype.(evlist{iEtype}) = cat(1,  X_ttype.(evlist{iEtype}), repmat(iEtype, [size(meanTFREv.(evlist{iEtype}),3), 1]));
                X_elec.(evlist{iEtype}) = cat(1, X_elec.(evlist{iEtype}), 10*iSub+iElec);
                X_subj.(evlist{iEtype}) = cat(1, X_subj.(evlist{iEtype}), iSub);
                
                avgLFP.(evlist{iEtype}) = cat(1, avgLFP.(evlist{iEtype}), meanEpochs);
                stdLFP.(evlist{iEtype}) = cat(1, stdLFP.(evlist{iEtype}), stdEpochs);
                
                [spectrum,ntaper,freqoi] = ft_specest_mtmfft(meanEpochs,[1:length(meanEpochs)]/srate,'freqoi',params.freqoiForAvgSpec,'taper','hanning');
                spectEv.(evlist{iEtype}) = cat(1, spectEv.(evlist{iEtype}), spectrum);
                nEvents.(evlist{iEtype}) = cat(1, nEvents.(evlist{iEtype}), nEpochs);
                chCnt(iEtype) = chCnt(iEtype) + 1;
                
                if iiCond == 3
                    disp(iiCond);
                end
                if isempty(X_subj_ch.(evlist{iEtype})(:))
                    warning('no data')
                    continue
                end
                GammaBehavTable_new = table(X_subj_ch.(evlist{iEtype})(:),X_elec_ch.(evlist{iEtype})(:),GammaBand_trial_all.(evlist{iEtype})(:),...
                    BehavD2D.(evlist{iEtype})(:),BehavD2G.(evlist{iEtype})(:),opt_path_impv12_trial.(evlist{iEtype})(:),...
                    'VariableNames',{'subj',  'elec', 'Gamma','distance2D','distance2G','opt_path_impv12'});
                GammaBehavTable.(evlist{iEtype}) = [GammaBehavTable.(evlist{iEtype}); GammaBehavTable_new];
                
            end
            
            
        end
        
        size(X_subj_ch.(evlist{iEtype})(:))
        size(X_elec_ch.(evlist{iEtype})(:))
        size(GammaBand_trial_all.(evlist{iEtype})(:))
        size(BehavD2D.(evlist{iEtype})(:))
        size(BehavD2G.(evlist{iEtype})(:))
        size(opt_path_impv12_trial.(evlist{iEtype})(:))
        
    end
    
end


% elec
clear data
filename = sprintf('CenteredAvgEvents_allElectrodesSubj_%s_%d.mat',roi,impvBlk2);
save(fullfile(savedir, 'population', filename),'freqoi','spectEv','nEvents','chCnt','avgLFP','stdLFP',...
    'meanTFREv','params','subjects','evlist','X_subj','X_elec','GammaBehavTable'); close all;


filename = sprintf('CenteredAvgEvents_allElectrodesSubj_%s_%d_TFRMAT.mat',roi,impvBlk2);
save(fullfile(savedir, 'population', filename),'freqoi','spectEv','nEvents','chCnt','avgLFP','stdLFP',...
    'meanTFREv','params','subjects','evlist','X_subj','X_elec','ALL_TFR_MAT','-v7.3'); close all;


end

% v2: fine-tuned twins for each contrasts (taking condition mean RTs into consideration)
% v3: updated to use sq_table instead of binfo
% v4: modified for testing interactions involving impv/worse 
% v5: stats within etype (shuffled on time-dimension); rt-outlier rm; 
% v6: grand-average - normalized RR, edited to focus analysis on improved
% paths

clc; clear all;  
%% define paths
addpath(genpath('C:\Users\mgeva\Documents\GitHub\external\eeglab2023.1'))
addpath(genpath('C:\Users\mgeva\Documents\GitHub\Kamin_iEEG_MAZE_codes'));
addpath(genpath('C:\Users\mgeva\Documents\GitHub\external\KK_useful\useful'));
addpath('C:\Users\mgeva\Documents\GitHub\external\fieldtrip-20230613\fieldtrip-20230613\external\brewermap')
maze_set_path; 

savedir = [dirs.ripplerate_ER 'traces_n7_v5/'];
%% ### Ideally, only edit this part ###
subjects = {'s01', 's06', 's08', 's15', 's16', 's17', 's18'}; 
roi = 'HPC' ;     
nperm = 250; 
grand_avg = 1; 
block_impv_grand_avg = 1; % only improved paths in grand average
block_effx = 1; 
block_impv_intx = 0; 
b3 = 0; 

%% color setting
% cmap = getPyPlot_cMap('Dark2');
cmap = brewermap(8,'BrBg');

plot_pars.col1  = cmap(5, :); 
plot_pars.col2  = cmap(8, :); 
plot_pars.ctrlcol = [0.6 0.6 0.6];

%% 
% events = {'N', 'D', 'O' 'G', 'X', 'E','mD', 'mO', 'mG'}; 
% events = {'mD', 'X', 'E', 'G'};  
% events = {'mD', 'G', 'E', 'X'};  
% events = { 'preG'};  
% events = {'G', 'E', 'X','mD', 'preG'};  
events = {'X','G', 'N','E', 'preG','mD'};  

for iEv = 1: length(events)  
    %% set condition names and twins
    eoi =  events{iEv};    

    switch eoi 
        case 'N' 
            trl_twin = [0.5 0.5]; 
        case 'D' 
            trl_twin = [0.5 0.5];            
        case 'O' 
            trl_twin = [0.35 0.7];  
        case 'G' 
            trl_twin = [0.5 0.5];  
        case 'X' 
            trl_twin = [0.2 1]; 
        case 'E'
            trl_twin = [0.8 0.2]; 
%         case 'mX' 
%             trl_twin = [1.3 0.2];            
        case 'mD' 
            %trl_twin = [0.85 0.2];            
            trl_twin = [0.7 0.2];            
%         case 'mO' 
%             trl_twin = [1.0 0.2]; 
        case 'mG' 
            trl_twin = [1.0 0.2];  
%         case 'preO' 
%             trl_twin = [0.6 0.1]; 
        case 'preG' 
            % trl_twin = [0.35 0.2];
            trl_twin = [0.5 0.5];            
     end
    
    %% start subj loop to aggregate data for traces
    % initialize agg structs
    RR_agg.swr = []; 
    RR_agg.block = []; 
    RR_agg.NMimpv12 = []; 
    RR_agg.NMimpv13 = []; 
    RR_agg.NMimpv123 = []; 
    RR_agg.elec = []; 
    RR_agg.subj = []; 
    RR_agg.baserr = []; 
    RR_agg.gdist = []; 
    RR_agg.ddist = []; 

    for iSub = 1: length(subjects) 
        subj = subjects{iSub};  
        roichans = get_roi_chans(subj, roi); 
       
         %% load data
        % Loading ripplerate 
        % RR data: for each square type, psth time-locked to each event occurance 
        % all events that in the sq_table 
        RR_eoi = load([dirs.ripplerate_ER, '/rr_', eoi, '/', subj '_RippleRate_ER_', eoi, '.mat']); % count per sec

        % behavioral 
        load([dirs.banal, 'sq_table/', subj '_sq_table.mat']); 
        [~, ~, rm_bool] = trim_mean(sq_table.rt , 3);        
        
        if strcmp(eoi(1), 'm')
            tmplist = strcmp(sq_table.square, eoi(end)); 
        elseif block_impv_grand_avg
            tmplist = strcmp(sq_table.square, eoi) & sq_table.opt_path_impv12 == 1; % GrandAverage
        else
            tmplist = strcmp(sq_table.square, eoi); % Orig code
        end 
        
        keeplist = tmplist; 
        keeplist = (tmplist & ~rm_bool); 
        keeplist = find(keeplist); 
        if strcmp(eoi(end), 'D') % D or mD
            [matchingD_ridx] =  func_reduce_D_v2 (subj); % reducte, match with RR
            rr_filter = ismember(keeplist, matchingD_ridx); 
            keeplist = keeplist(rr_filter);
        end
        
        rr_ridx = find(ismember(find(tmplist),keeplist));
        eoi_ridx = keeplist;
        if length(rr_ridx) ~= length(eoi_ridx)
            error('check number of row indices');
        end
        
        %% RR tpnts: defining tpnts & edges within a trial
        % RR are obtained for [1.5 1.5] window and this step allows for 
        % analayses limited to the twin 
        RR_tpnts =RR_eoi.edges(1,:) - mean(RR_eoi.edges(1,1:2)) - RR_eoi.twin(1);
        RR_tpnts = RR_tpnts(2: end); 
        RR_edges =RR_eoi.edges(1,:) - RR_eoi.edges(1, 1)- RR_eoi.twin(1); 

        tmarks = [trl_twin(1)*(-1) 0 trl_twin(2)];
        RR_on =max(find(RR_tpnts'<tmarks(1)));
        RR_off = min(find(RR_tpnts' > tmarks(end)));
        RR_tpnts = RR_tpnts(RR_on: RR_off); 
        %% pick roi channels & rectify names if needed 
        subjROI = roichans;
        for i = 1: length(subjROI)
            subjROI{i}(strfind(subjROI{i}, ' '))=[];
            subjROI{i}(strfind(subjROI{i}, '-'))=[];
        end        
        chans = fieldnames(RR_eoi.EventRR); 
        cidx = find(ismember(chans, subjROI));

        %% get task timing to use for baseline task rr
        load([dirs.ripple, '/', subj, '/' subj '_ripples.mat']); 
        task_onset = sq_table.time(1); 
        task_offset = sq_table.time(end)+1; 
        
        %% loop over subject electrodes
        for iChan = 1: length(cidx)
            % baseline ripple rate
            ripple_peaks = ripples.(chans{cidx(iChan)}).peak; 
            ripple_peaks = ripple_peaks(ripple_peaks>task_onset & ripple_peaks <task_offset);          
            baserate = length(ripple_peaks) / (task_offset - task_onset); 
            
            RR_traces_eoi = RR_eoi.EventRR.(chans{cidx(iChan)}); 
            %             if sum(strcmp(sq_table.square, eoi(end))) ~= size(RR_eoi.EventRR.(chans{cidx(iChan)}), 1) 
            if sum(strcmp(sq_table.square, eoi)) ~= size(RR_eoi.EventRR.(chans{cidx(iChan)}), 1)
                if (strcmp(eoi,'mD') || strcmp(eoi,'mG') )&& ~sum(strcmp(sq_table.square, eoi(end))) ~= size(RR_eoi.EventRR.(chans{cidx(iChan)}), 1)
                    disp(eoi)
                else
                    error('check number of event types in sq_table');
                end
            end
            
            %% aggregate over channels & subjs 
            RR_agg.swr      = cat(1, RR_agg.swr, RR_traces_eoi(rr_ridx, :)); 
            RR_agg.block    = cat(1, RR_agg.block, sq_table.block(eoi_ridx)); 
            RR_agg.NMimpv12 = cat(1, RR_agg.NMimpv12, sq_table.opt_path_impv12(eoi_ridx)); 
            RR_agg.NMimpv13 = cat(1, RR_agg.NMimpv13, sq_table.opt_path_impv13(eoi_ridx)); 
            if ismember(3, sq_table.block(eoi_ridx))
              RR_agg.NMimpv123 = cat(1, RR_agg.NMimpv123, sq_table.opt_path_impv123(eoi_ridx));       
            else
              RR_agg.NMimpv123 = cat(1, RR_agg.NMimpv123, sq_table.opt_path_impv12(eoi_ridx)); 
            end             
            RR_agg.baserr     = cat(1, RR_agg.baserr, repmat(baserate, [length(eoi_ridx) 1])); 
            RR_agg.elec     = cat(1, RR_agg.elec, repmat([iSub*10+iChan], [length(eoi_ridx) 1])); 
            RR_agg.subj     = cat(1, RR_agg.subj, repmat([iSub], [length(eoi_ridx) 1])); 
            RR_agg.gdist    = cat(1, RR_agg.gdist, sq_table.d2goal(eoi_ridx)); 
            RR_agg.ddist    = cat(1, RR_agg.ddist, sq_table.d2dec(eoi_ridx)); 
        end
        clear RR_eoi eoi_ridx 
    end           
    %% reduce RR_agg to d2goal > 2 , d2dec > 2 and to event-specific time window
    vars = fieldnames(RR_agg); 
    vars =  setdiff(vars, {'gdist', 'ddist'}); 
    if strcmp(eoi(end), 'D') 
        for ii = 1:length(vars)
            RR_agg.(vars{ii}) = RR_agg.(vars{ii}) (abs(RR_agg.gdist) > 2, :); 
        end         
    elseif strcmp(eoi(end), 'O') 
        for ii = 1:length(vars)
            RR_agg.(vars{ii}) = RR_agg.(vars{ii}) (abs(RR_agg.gdist) > 2 & RR_agg.ddist ~=1 &  RR_agg.ddist ~=2, :); 
        end
    elseif strcmp(eoi(end), 'G') 
        for ii = 1:length(vars)
            RR_agg.(vars{ii}) = RR_agg.(vars{ii}) (abs(RR_agg.ddist) > 2, :); 
        end
    end 
    tracedata = RR_agg.swr(:, RR_on: RR_off); 
    %% Trace Grand Average 
    if grand_avg
        if b3
            savename = ['RR_grandavg_' eoi '\'];
        elseif block_impv_grand_avg
            savename = ['RR_grandavg_impvOnly_b12_' eoi ];
        else
            savename = ['RR_grandavg_' eoi '_b12_2tail'];
        end
        
        %% statistical testing for grand average
        % create permutation distribution [shuffle on time dimension]
        permRR = nan(nperm, size(tracedata, 2));     
        permtrace = nan([size(tracedata), nperm]);
        for iPerm = 1: nperm
            randcut = randsample([round(size(tracedata, 2)*0.15) : round(size(tracedata, 2)*0.85)], size(tracedata, 1), true);            
            for ii = 1: size(tracedata, 1)
                permtrace(ii, :, iPerm) = cat(2, tracedata(ii, randcut(ii)+1:end), tracedata(ii, 1: randcut(ii))); 
            end 
            permRR(iPerm, :) = mean(permtrace(:, :, iPerm), 1);         
        end
        
        ztrace = (mean(tracedata, 1) - mean(permRR, 1)) ./ std(permRR, 0, 1);
        ztracethresh = ztrace; 
        ztracethresh(abs(ztracethresh) < norminv(1-0.025)) = 0 % 2-tailed
        % apply cluster thresholding 
        clustinfo = bwconncomp(ztracethresh);
        clust_info = cellfun(@numel,clustinfo.PixelIdxList);  
        perm_clust_info = []; 
        for iPerm = 1: nperm
            permz = (permRR(iPerm, :) - mean(permRR, 1)) ./ std(permRR, 0, 1);        
            perm_tracethres = abs(permz) > norminv(1-0.025);
            permclustinfo = bwconncomp(perm_tracethres);
            perm_clust_info = [perm_clust_info cellfun(@numel,permclustinfo.PixelIdxList)];  
        end
        prctile(perm_clust_info, 95);
        whichclusters2remove = find(clust_info<prctile(perm_clust_info, 95));
        % remove clusters
        for i=1:length(whichclusters2remove)
            ztracethresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
        end
        sig_masked = double(ztracethresh ~=0);

        %% plot: grand average 
        ylimunits = [0:0.5:1];
        figure('position', [400, 400, 430, 480],'Color','w');

        tracedata = RR_agg.swr(:, RR_on: RR_off); 
        % traces with std
        h = plot(RR_tpnts, mean(tracedata),'Color', plot_pars.col1, 'LineWidth', 2);  hold on; 
        jbfill(RR_tpnts, mean(tracedata)+std(tracedata, 0, 1)/sqrt(size(tracedata,1)-1), ...
            mean(tracedata)-std(tracedata, 0, 1)/sqrt(size(tracedata,1)-1), plot_pars.col1, plot_pars.col1,  1, 0.2); 
%         ylims = ylim; 
%         a = dsearchn(ylimunits', ylims(1));     b = dsearchn(ylimunits', ylims(2));
%         if ylimunits(a) > ylims(1), ylims(1) = ylimunits(a-1); else, ylims(1) = ylimunits(a); end 
%         if ylimunits(b) < ylims(2), ylims(2) = ylimunits(b+1); else, ylims(2) = ylimunits(b); end     
        ylim([0 1]); ylims = ylim; 
        if ~isempty(find(sig_masked == 1))
            clear clustinfo 
            clustinfo = bwconncomp(sig_masked);
            for iSigBump = 1: clustinfo.NumObjects
                BumpCntr = clustinfo.PixelIdxList{iSigBump}';
                jbfill(RR_tpnts(BumpCntr), repmat([ylims(2)], [1 length(RR_tpnts(BumpCntr))]),  repmat([ylims(2)*0.98], [1 length(RR_tpnts(BumpCntr))]),[170 21 27]/255, [170 21 27]/255,  1); 
            end 
        end 
        % plot look setting
        xticks(tmarks); 
        xticklabels(tmarks); 
        yticks([ylims(1) : 0.5: ylims(2)]); 
        xlim([tmarks(1) tmarks(end)]); 
        xlabel('time (sec)')
        ylabel('ripple rate (count/s)')   
        title(['grand average ',eoi])
        set(gca, 'FontSize', 16); 
        line ([0 0], ylims, 'LineWidth', 1.5, 'Color', [0.4 0.4 0.4], 'LineStyle', '--');
        % save figure 
        saveas(gcf, [savedir, savename], 'png');        
    end
    %% Set up data for GLMM
    Y = tracedata;
    Xtable = table(RR_agg.block, RR_agg.NMimpv12, RR_agg.NMimpv13,RR_agg.NMimpv123, RR_agg.elec, RR_agg.subj, 'VariableNames', {'block', 'impv12', 'impv13', 'impv123','elec', 'subj'});     
    %%  Block effect
    if block_effx
        if b3 == 0
            Xtable(RR_agg.block == 3,:) = []; 
            Y(RR_agg.block == 3,:) =[]; 
            savename = ['RR_block_effx_b12_' eoi ]; 
            vars = {'block'}; 
            permvar = {'block'};            
            rr_data.c1 = Y(Xtable.block == 1, :);% target condition is coded as 2
            rr_data.c2 = Y(Xtable.block == 2, :); % Y(Xtable.ttype == 1, :);
        else
            Xtable.repetition = (Xtable.block > 1); % repetition vs not
            savename = ['RR_block_effx_b123_' eoi ]; 
            vars = {'repetition'}; 
            permvar = {'repetition'};   
            rr_data.c1 = Y(Xtable.repetition == 0, :);% target condition is coded as 2
            rr_data.c2 = Y(Xtable.repetition == 1, :); % Y(Xtable.ttype == 1, :);
        end 
        
        %% statistical testing for b1-b2 difference
        %% run GLMM
        % [real_stat, perm_stat, sig_masked, sig_masked_noclust] = trace_GLMM_wPerm(Y, Xtable, nperm, vars, permvar, 0); % trace_GLMM_wPerm(Y, Xtable, nperm, permvars, intx)
        switch eoi
            case {'mD', 'preG'}
                [real_stat, perm_stat, sig_masked, sig_masked_noclust] = trace_GLMM_wPerm(Y(:, RR_tpnts<0), Xtable, nperm, vars, permvar, 0); 
            case {'G', 'E', 'X'}
                [real_stat, perm_stat, sig_masked, sig_masked_noclust] = trace_GLMM_wPerm(Y(:, RR_tpnts>0), Xtable, nperm, vars, permvar, 0); 
        end 

        zstat = (real_stat - mean(perm_stat, 2)) ./ std (perm_stat, 0, 2);
        save([savedir, savename, 'zstat.mat'], 'zstat','vars','Xtable','RR_tpnts','nperm')
        
        sig_masked(isnan(sig_masked)) = 0;
        
        rr_parc.edges = RR_edges(dsearchn(RR_edges', tmarks(1)'): dsearchn(RR_edges', tmarks(end)')); 
        rr_pars.tpnts = RR_tpnts;
        rr_pars.sigmarks = nan(length(sig_masked)); % same length as Y in time dim;
        rr_pars.yrange = [0 1];
        plot_pars.col2 = plot_pars.col1*0.6;   % set color for 2nd block
        plot_trace (plot_pars, rr_data, rr_pars, tmarks); 
        yticks([rr_pars.yrange(1) : 0.5: rr_pars.yrange(2)]); 
        
        if ~isempty(find(sig_masked == 1))
            clear clustinfo 
            clustinfo = bwconncomp(sig_masked);
            for iSigBump = 1: clustinfo.NumObjects
                BumpCntr = clustinfo.PixelIdxList{iSigBump}';
                jbfill(RR_tpnts(BumpCntr), repmat([rr_pars.yrange(2)], [1 length(RR_tpnts(BumpCntr))]),  repmat([rr_pars.yrange(2)*0.95], [1 length(RR_tpnts(BumpCntr))]),[170 21 27]/255, [170 21 27]/255,  0.1); 
            end 
        end
        title(['block_effx ',eoi])

        saveas(gcf, [savedir, savename], 'png');
        
        args = {[1:length(zstat)], [1 2], [zstat, zstat]'}; %  
        map_handle = pcolor(args{:}); colormap(parula)
        set(map_handle, 'EdgeColor', 'none'); 
        shading interp;           
        set(gca,'Box', 'off'); hold on;     
        colorbar ('YTick', [-1.96 -1.65 0 1.65 1.96]); caxis([-1.95 1.96]);        
        saveas(gcf, [savedir, savename, 'zbar'], 'png');      

    end 
    
    if block_impv_intx
        if b3 == 0
            if strcmp(eoi,'preG')
                rmvRow = find(Xtable.block == 3);
                Xtable(rmvRow,:) = [];
                Y(rmvRow,:) = [];
            else
                Xtable(RR_agg.block == 3,:) = [];
                Y(RR_agg.block == 3,:) =[];
            end
            vars = {'block', 'impv12'}; 
            permvar = {'block', 'impv12'}; 
            savename = ['RR_block_by_impv_intx_' eoi ];
        else
            Xtable.repetition = (Xtable.block > 1); % repetition vs not
            vars = {'repetition', 'impv123'}; 
            permvar = {'repetition', 'impv123'}; 
            savename = ['RR_block_by_impv_intx_' eoi '_b123'];
        end
        [real_stat, perm_stat, sig_masked, sig_masked_noclust] = trace_GLMM_wPerm(Y, Xtable, nperm, vars, permvar, 1); % trace_GLMM_wPerm(Y, Xtable, nperm, permvars, intx)
        zstat = (real_stat - mean(perm_stat, 2)) ./ std (perm_stat, 0, 2);
        zstat_intx = zstat;
        sig_masked_noclust_intx = sig_masked_noclust;
        
        vars = {'block'}; 
        permvar = {'block'};
        rmvRow = find(Xtable.impv12 == 0);
        Y1 = Y;
        Xtable1 = Xtable;
        Xtable1(rmvRow,:) = [];
        Y1(rmvRow,:) =[]; 
        [real_stat12, perm_stat12, sig_masked12, sig_masked_noclust12] = trace_GLMM_wPerm(Y, Xtable, nperm, vars, permvar, 0); % trace_GLMM_wPerm(Y, Xtable, nperm, permvars, intx)
        zstat12 = (real_stat - mean(perm_stat, 2)) ./ std (perm_stat, 0, 2);

        newA4figure("")
        for ii_a = 1:3
            
            if ii_a == 1
                axes('position',[0.1 0.5 0.3 0.2])
                if b3 == 0
                    rr_data.c1 = Y(Xtable.block == 1 & Xtable.impv12 == 1, :);% target condition is coded as 2
                    rr_data.c2 = Y(Xtable.block == 2 & Xtable.impv12 == 1, :); % Y(Xtable.ttype == 1, :);
                else
                    rr_data.c1 = Y(Xtable.repetition == 0 & Xtable.impv123 == 1, :);% target condition is coded as 2
                    rr_data.c2 = Y(Xtable.repetition == 1 & Xtable.impv123 == 1, :); % Y(Xtable.ttype == 1, :);
                end
            elseif ii_a == 2
                clear rr_data
                axes('position',[0.1 0.1 0.3 0.2])
                if b3 == 0
                    rr_data.c1 = Y(Xtable.block == 1 & Xtable.impv12 == 0, :);% target condition is coded as 2
                    rr_data.c2 = Y(Xtable.block == 2 & Xtable.impv12 == 0, :); % Y(Xtable.ttype == 1, :);
                else
                    rr_data.c1 = Y(Xtable.repetition == 0 & Xtable.impv123 == 0, :);% target condition is coded as 2
                    rr_data.c2 = Y(Xtable.repetition == 1 & Xtable.impv123 == 0, :); % Y(Xtable.ttype == 1, :);
                end
                
            elseif ii_a == 3
                % add another plot to specifically compare block 1-2 for improved
                % path
                clear rr_data
                axes('position',[0.5 0.5 0.3 0.2])
                if b3 == 0
                    rr_data.c1 = Y(Xtable.block == 1 & Xtable.impv12 == 1, :);% target condition is coded as 2
                    rr_data.c2 = Y(Xtable.block == 2 & Xtable.impv12 == 1, :); % Y(Xtable.ttype == 1, :);
                else
                    rr_data.c1 = Y(Xtable.repetition == 0 & Xtable.impv123 == 1, :);% target condition is coded as 2
                    rr_data.c2 = Y(Xtable.repetition == 1 & Xtable.impv123 == 1, :); % Y(Xtable.ttype == 1, :);
                end
                
                sig_masked_noclust = sig_masked_noclust12;
                zstat = zstat12;
                
            end
            
            
            
            rr_parc.edges = RR_edges(dsearchn(RR_edges', tmarks(1)'): dsearchn(RR_edges', tmarks(end)'));
            rr_pars.tpnts = RR_tpnts;
            rr_pars.sigmarks = sig_masked_noclust; % same length as Y in time dim;
            rr_pars.yrange = [0 1];
            
            % traces with std
            plot(rr_pars.tpnts, mean(rr_data.c1),'Color', plot_pars.col1, 'LineWidth', 2);  hold on;
            plot(rr_pars.tpnts, mean(rr_data.c2),'Color', plot_pars.col1*0.6, 'LineWidth', 2);  hold on;
            jbfill(rr_pars.tpnts, mean(rr_data.c1)+std(rr_data.c1, 0, 1)/sqrt(size(rr_data.c1,1)-1), ...
                mean(rr_data.c1)-std(rr_data.c1, 0, 1)/sqrt(size(rr_data.c1,1)-1), plot_pars.col1, plot_pars.col1,  1, 0.2);
            jbfill(rr_pars.tpnts, mean(rr_data.c2)+std(rr_data.c2, 0, 1)/sqrt(size(rr_data.c2,1)-1), ...
                mean(rr_data.c2)-std(rr_data.c2, 0, 1)/sqrt(size(rr_data.c2,1)-1), plot_pars.col1*0.6,  plot_pars.col1*0.6,  1, 0.2);  hold on;
            
            % sigmark shades
            sig_masked = rr_pars.sigmarks;
            sig_clusts = find(~isnan(sig_masked))';
            sig_masked(isnan(sig_masked)) = 0;
            clustinfo = bwconncomp(sig_masked);
            for iSigBump = 1: clustinfo.NumObjects
                BumpCntr = clustinfo.PixelIdxList{iSigBump}';
                jbfill(rr_pars.tpnts(BumpCntr), repmat([15], [1 length(rr_pars.tpnts(BumpCntr))]),  repmat([-1], [1 length(rr_pars.tpnts(BumpCntr))]),[253 210 110]/255, [253 210 110]/255,  1, 0.2);
            end
            line ([0 0], rr_pars.yrange, 'LineWidth', 1.5, 'Color', [0.4 0.4 0.4], 'LineStyle', '--');
            
            % plot look setting
            xticks(tmarks);
            xticklabels(tmarks);
            xlim([tmarks(1) tmarks(end)]);
            ylim(rr_pars.yrange);
            xlabel('time (sec)')
            ylabel('ripple rate (count/s)')
            set(gca, 'FontSize', 16)
            
            if ii_a == 2; continue;
            if ii_a == 1
                axes('position',[0.1 0.7 0.3 0.05])
            elseif ii_a == 3
                axes('position',[0.5 0.7 0.3 0.05])
            end
                args = {[1:length(zstat)], [1 2], [zstat, zstat]'}; %
                map_handle = pcolor(args{:}); colormap(parula)
                set(map_handle, 'EdgeColor', 'none');
                shading interp;
                set(gca,'Box', 'off'); hold on;
                set(gca, 'xtick',[],'ytick',[])
                axes('position',[0.35 0.65 0.1 0.1])
                colorbar ('YTick', [-1.96 -1.65 0 1.65 1.96]); caxis([-1.95 1.96]);
                axis off
            end
            
        end
        
        saveas(gcf, [savedir, savename], 'png');  

    end 
end 

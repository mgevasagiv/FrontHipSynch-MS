% This is a master code for analyzing behavioral data from MAZE task 
% Meant to be run by sections as needed 
% Input:log.csv 
% OUTPUT (saves): 
% 1. logtable (condensed with relevant info only) - does this with manual
% edition for subjects whose logs did not have maze info
% 2. generates and saves a maze_difficulty table 
% 3. runs func_nave_impv
% 4. group-level summary plots: navigation performance over b1 & b2 
% 5. group-level summary plots by maze difficulty 
% 6. get number of matching dp points for impv and non-impv mazes -- each
% subject

%% set dir & subjs
maze_set_path; 
logdir = [dirs.banal 'logtables/']; 
subjs = {'s01', 's06', 's08', 's09', 's12', 's14', 's15', 's16', 's17', 's18'}; 
% subjs = {'s18'}; 
%% 1. save logtables 
for iSub = 1: length(subjs)
    subj = subjs{iSub}; 
    if str2num(subj(strfind(subj, 's')+1:end)) < 10, beh_type = 1; else, beh_type = 2; end 
    logtable = get_logtable (subj, beh_type, dirs.behavioral); 
    save(fullfile(logdir, [subj, '_logtable.mat']), 'logtable'); 
end 
%% MAZE names (world)were logged in raw data from s08 on.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NOT IN IR SUBJECTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Therefore, maze names were manually added for s01, s06, s16, s17; saved as *new.xslx[Radhika] 
% Here, just converting/saving them into .mat for consistency 
editsubjs = {'s01', 's06'}; %, 's16', 's17'};
for iSub = 1: length(editsubjs)
    subj = editsubjs{iSub}; 
    logtable = readtable(fullfile(logdir, [subj '_new.xlsx']));
    save(fullfile(logdir, [subj, '_logtable.mat']), 'logtable'); 
    clear logtable; 
end 
%% 2. Create & save maze_difficulty.mat -- leaving it here for the record on how that was done
% tmp_difftable = readtable('/Users/kaminkim/Documents/projects/iEEG_MAZE/analysis/behav/difficulty/difficulty.xlsx');
% tmp_worldtable = cell2table(unique(logtable.world));
% maze_difficulty = [tmp_worldtable tmp_difftable]; 
% save(fullfile('/Users/kaminkim/Documents/codes/iEEG_MAZE/behav', 'maze_difficulty.mat'), 'maze_difficulty'); 
%% 3. Classifying mazes by performance & difficulty 
outdir =  [dirs.banal 'nav_impv_v3/']; 
for iSub = 1:length(subjs)
    subj = subjs{iSub}; 
    func_nav_impv_v3 (subj, logdir, outdir);  % saves a table with 30 mazes - trial idx & navigatin record in b1 & b2
    % need to additionally save 'which decision points are comparable
    % across blocks (i.e., the same decision points)' 
end 

%% MAZE0_get_beh_indv_v2 should be run here (after get_logtable, func_nav_impv, and MAZE1_get_trigger)
% problem with this is that some 'D's that are counted in func_nav_impv_v3
% get 'nulled' by MAZE_bet_beh_indv_v2 --> which makes func_reduce_v2 break
% (s09) 

%% 4. Group summary of behavioral changes in B1 & B2
savedir = dirs.banal;
for iSub = 1: length(subjs)
    subj = subjs{iSub}; 
    load (fullfile([dirs.banal 'nav_impv_v3/'], [subj '_maze_nav_impv.mat']));
    bdata.subj{iSub,1} = subj;   
    bdata.nmove_b1(iSub)  = mean(maze_nav_impv.nMoves_b1); 
    bdata.nmove_b2(iSub)  = mean(maze_nav_impv.nMoves_b2);
    bdata.t2complete_b1(iSub)  = mean(maze_nav_impv.TotalTime_b1);
    bdata.t2complete_b2(iSub)  = mean(maze_nav_impv.TotalTime_b2);
    bdata.t2dec_b1(iSub)  = nanmean(maze_nav_impv.DecisionTime_b1);  
    bdata.t2dec_b2(iSub)  = nanmean(maze_nav_impv.DecisionTime_b2);
    bdata.optsteps_b1(iSub)  = nanmean(maze_nav_impv.OptimalStepOffset_b1);  
    bdata.optsteps_b2(iSub)  = nanmean(maze_nav_impv.OptimalStepOffset_b2);
    bdata.optpath_b1(iSub)  = nanmean(maze_nav_impv.OptimaPathOffset_b1);  
    bdata.optpath_b2(iSub)  = nanmean(maze_nav_impv.OptimaPathOffset_b2);

end 

% savename = ['bsummary_n' num2str(length(subjs))]; 
% table2save = table(bdata.subj(:), bdata.nmove_b1(:), bdata.nmove_b2(:), ...
%     bdata.t2complete_b1(:), bdata.t2complete_b2(:), bdata.t2dec_b1(:), bdata.t2dec_b2(:), ...
%     'VariableNames', {'subj', 'nmove_b1', 'nmove_b2', 't2c_b1', 't2c_b2', 't2d_b1', 't2d_b2'}); 
% writetable(table2save, fullfile(dirs.banal, [savename '.csv']));

savename = ['NavImprovement_n' num2str(length(subjs))]; 
figure('units','normalized','outerposition',[0.2 0.2 0.6 0.4],'Color','w');
subplot(1, 5, 1);bplot_overblocks(bdata.nmove_b1(:), bdata.nmove_b2(:), 'number of moves', 'count');
subplot(1, 5, 2);bplot_overblocks(bdata.t2complete_b1(:)/1000, bdata.t2complete_b2(:)/1000, 'time to complete', 'sec');
subplot(1, 5, 3);bplot_overblocks(bdata.t2dec_b1(:), bdata.t2dec_b2(:),'decision time', 'sec');
subplot(1, 5, 4);bplot_overblocks(bdata.optsteps_b1(:), bdata.optsteps_b2(:), 'offset from optimal #steps', 'count');
subplot(1, 5, 5);bplot_overblocks(bdata.optpath_b1(:), bdata.optpath_b2(:),'offset from optimal paths', 'count');
legend(bdata.subj); legend boxoff; 
saveas(gcf,fullfile(savedir, [savename '.png']),'png'); 

%% 5. Group summary by maze difficulty
diffs = {'easy', 'medium', 'hard'};
for iSub = 1: length(subjs)
    subj = subjs{iSub}; 
    load (fullfile([dirs.banal 'nav_impv_v3/'], [subj '_maze_nav_impv.mat']));
    for iDiff = 1: length(diffs)
        midx = find(strcmp(maze_nav_impv.difficulty, diffs{iDiff})); 
        NavPerfs.TotalTime.(diffs{iDiff})(iSub, :) = nanmean([maze_nav_impv.TotalTime_b1(midx) maze_nav_impv.TotalTime_b2(midx)]);
        NavPerfs.nMoves.(diffs{iDiff})(iSub, :) = nanmean([maze_nav_impv.nMoves_b1(midx) maze_nav_impv.nMoves_b2(midx)]); 
        NavPerfs.DecisionTime.(diffs{iDiff})(iSub, :) = nanmean([maze_nav_impv.DecisionTime_b1(midx) maze_nav_impv.DecisionTime_b2(midx)]); 
    end 
end 

savename = ['DecitionTime_by_Difficulty_n' num2str(length(subjs))]; 
figure('units','normalized','outerposition',[0.2 0.2 0.6 0.4],'Color','w');
subplot(1, 3, 1);bplot_overblocks(NavPerfs.DecisionTime.easy(:, 1)/1000, NavPerfs.DecisionTime.easy(:, 2)/1000, 'easy', 'sec');
subplot(1, 3, 2);bplot_overblocks(NavPerfs.DecisionTime.medium(:, 1)/1000, NavPerfs.DecisionTime.medium(:, 2)/1000, 'medium', 'sec');
subplot(1, 3, 3);bplot_overblocks(NavPerfs.DecisionTime.hard(:, 1)/1000, NavPerfs.DecisionTime.hard(:, 2)/1000, 'hard', 'sec');
legend(bdata.subj); legend boxoff; 
saveas(gcf,fullfile(savedir, [savename '.png']),'png'); 

savename = ['TotalTime_by_Difficulty_n' num2str(length(subjs))]; 
figure('units','normalized','outerposition',[0.2 0.2 0.6 0.4],'Color','w');
subplot(1, 3, 1);bplot_overblocks(NavPerfs.TotalTime.easy(:, 1)/1000, NavPerfs.TotalTime.easy(:, 2)/1000, 'easy', 'sec');
subplot(1, 3, 2);bplot_overblocks(NavPerfs.TotalTime.medium(:, 1)/1000, NavPerfs.TotalTime.medium(:, 2)/1000, 'medium', 'sec');
subplot(1, 3, 3);bplot_overblocks(NavPerfs.TotalTime.hard(:, 1)/1000, NavPerfs.TotalTime.hard(:, 2)/1000, 'hard', 'sec');
legend(bdata.subj); legend boxoff; 
saveas(gcf,fullfile(savedir, [savename '.png']),'png'); 

savename = ['nMove_by_Difficulty_n' num2str(length(subjs))]; 
figure('units','normalized','outerposition',[0.2 0.2 0.6 0.4],'Color','w');
subplot(1, 3, 1);bplot_overblocks(NavPerfs.nMoves.easy(:, 1), NavPerfs.nMoves.easy(:, 2), 'easy', 'count');
subplot(1, 3, 2);bplot_overblocks(NavPerfs.nMoves.medium(:, 1), NavPerfs.nMoves.medium(:, 2), 'medium', 'count');
subplot(1, 3, 3);bplot_overblocks(NavPerfs.nMoves.hard(:, 1), NavPerfs.nMoves.hard(:, 2), 'hard', 'count');
legend(bdata.subj); legend boxoff; 
saveas(gcf,fullfile(savedir, [savename '.png']),'png'); 

%% 6. get number of matching dp points for impv and non-impv mazes -- each subject
nMatch = nan(length(subjs), 5); % all, DTimpv, DTnoimpv, NMimpv, NMnoimpv
for iSub = 1: length(subjs)
    subj = subjs{iSub}; 
    load([dirs.behavioral, '/', subj '_binfo.mat']); 
    load([dirs.banal, 'nav_impv_v2/', subj '_matching_decpnts.mat']);
    
    % get number of matching dp in each subject (from dpmatch & func_reducedD)
    % care about this to be aware of each subject's contribution to total
    % nTrials that go into analyses 
    [reducedD, rmD] = func_reduce_D (square.D, dpmatch); 
    nMatch(iSub, 1) =  length(reducedD)/2; 
      
    % classify by reduced dTime  
    DT_dpmatch = []; 
    for iMaze = 1: length(dpmatch)
        dEvt_idx_b1 = find(([square.D.block] == 1) & ([square.D.tidx] == dpmatch(iMaze).tidx_b1));
        dEvt_idx_b2 = find(([square.D.block] == 2) & ([square.D.tidx] == dpmatch(iMaze).tidx_b2));
        DT_dpmatch = cat(1, DT_dpmatch, ...
            cat(2, [square.D(dEvt_idx_b1(dpmatch(iMaze).matchingdp_b1 ==1)).t2dec]', ...
            [square.D(dEvt_idx_b2(dpmatch(iMaze).matchingdp_b2 ==1)).t2dec]'));
        clear dEvt_idx*
    end 
    nMatch(iSub, 2) = sum(DT_dpmatch(:, 1) - DT_dpmatch(:, 2) > 0); 
    nMatch(iSub, 3) = sum(DT_dpmatch(:, 1) - DT_dpmatch(:, 2) <= 0);
    figure; plot(DT_dpmatch', 'ko-', 'LineWidth', 1); xlim([0.5 2.5]); xticks([1 2]);  
    
    % classify by reduced nMoves 
    load (fullfile([dirs.banal 'nav_impv_v2/'], [subj '_maze_nav_impv.mat']));    
    NM_dpmatch = nan(length(dpmatch), 2); 
    for iMaze = 1: length(dpmatch)        
        n_mdpnts = sum(dpmatch(iMaze).matchingdp_b1 ==1); % dpnts across b1&b2 in this maze
        n_mdpnts = [n_mdpnts 0]; 
        if maze_nav_impv.nMoves_impv(strcmp(maze_nav_impv.world, dpmatch(iMaze).mazename)) == 1  
            NM_dpmatch(iMaze, :) = n_mdpnts; 
        else
            NM_dpmatch(iMaze, :) = fliplr(n_mdpnts); 
        end
    end    
    nMatch(iSub, 4:5) = sum(NM_dpmatch, 1); 
end 

%% nMoves  
figure
for iSub = 1: length(subjs)
    subj = subjs{iSub}; 
    load (fullfile([dirs.banal 'nav_impv_v3/'], [subj '_maze_nav_impv.mat']));    
    subplot(5,2,iSub)
    plot([maze_nav_impv.nMoves_b1 maze_nav_impv.nMoves_b2]', 'o-', 'LineWidth', 1,...
        'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.3 0.3 0.3],  'MarkerEdgeColor', [0.3 0.3 0.3]); 
    xlim([0.5 2.5]); xticks([1 2]);  
    xticklabels({'block1', 'block2'});
    ylim([0 50]); 
    set(gca, 'FontSize', 16); box off
%     saveas(gcf,fullfile(dirs.banal, ['subj_nMoves_' subj '.tiff']),'tiff'); 
end 

%% nMoves - group, impv vs worse, for paper
cmap = colormap(jet); 
cmap = cmap(round(linspace(1, 64, 10)), :);
agg_impv = [];
agg_same = []; 
agg_worse = [];
agg_diff_steps= [];
ratio = nan(length(subjs), 3); 

for iSub = 1: length(subjs)
    subj = subjs{iSub}; 
    load (fullfile([dirs.banal 'nav_impv_v3/'], [subj '_maze_nav_impv.mat']));  
    
    maze_nav_impv(isnan(maze_nav_impv.DecisionTime_b1),:) = []; 
    
    agg_impv = cat(1, agg_impv, [maze_nav_impv.nMoves_b1(maze_nav_impv.nMoves_b1 > maze_nav_impv.nMoves_b2), maze_nav_impv.nMoves_b2(maze_nav_impv.nMoves_b1 > maze_nav_impv.nMoves_b2), repmat(iSub, [sum(maze_nav_impv.nMoves_b1 > maze_nav_impv.nMoves_b2) 1])]); 
    agg_same = cat(1, agg_same, [maze_nav_impv.nMoves_b1(maze_nav_impv.nMoves_b1 == maze_nav_impv.nMoves_b2), maze_nav_impv.nMoves_b2(maze_nav_impv.nMoves_b1 == maze_nav_impv.nMoves_b2), repmat(iSub, [sum(maze_nav_impv.nMoves_b1 == maze_nav_impv.nMoves_b2) 1])]); 
    agg_worse = cat(1, agg_worse, [maze_nav_impv.nMoves_b1(maze_nav_impv.nMoves_b1 < maze_nav_impv.nMoves_b2), maze_nav_impv.nMoves_b2(maze_nav_impv.nMoves_b1 < maze_nav_impv.nMoves_b2), repmat(iSub, [sum(maze_nav_impv.nMoves_b1 < maze_nav_impv.nMoves_b2) 1])]); 
    
    agg_diff_steps = cat(1, agg_diff_steps, [maze_nav_impv.nMoves_b2-maze_nav_impv.nMoves_b1]); 

    ratio(iSub, :) = [sum(maze_nav_impv.nMoves_b1 > maze_nav_impv.nMoves_b2)/height(maze_nav_impv), ...
                    sum(maze_nav_impv.nMoves_b1 == maze_nav_impv.nMoves_b2)/height(maze_nav_impv), ...
                    sum(maze_nav_impv.nMoves_b1 < maze_nav_impv.nMoves_b2)/height(maze_nav_impv)] ;        
end  

newA4figure('maze_behav_ratios')
cmap_b = brewermap(3,'Blues');
for iSub = 1: length(subjs)
    subplot(2,5,iSub)
    h = pie(ratio(iSub, :),[0,0,1]); 
    
    if ~sum(ismember(ratio(iSub, :),0))
        h(1).FaceColor = cmap_b(2,:);  % blue
        h(3).FaceColor = cmap_b(1,:);  % blue
        h(5).FaceColor = [0.8 0.8 0.8];  % gray
    else
        h(1).FaceColor = cmap_b(2,:);  % Red
        h(3).FaceColor = [0.8 0.8 0.8];  % gray
    end
    
    title(sprintf('s%d',iSub))
%     for k = 1:length(h)
%      if isgraphics(h(k), 'text')
%         delete(h(k));
%      end
%     end
end

savename = ['MGS_ratioPerSubj_n', num2str(length(subjs))];
saveas(gcf,fullfile(dirs.banal, [savename '.tiff']),'tiff'); 

% cmap = colormap(jet); 
% cmap = cmap(round(linspace(1, 256, 10)), :);

% cmap = [0.55294118, 0.82745098, 0.78039216;
%  1.,         1.,         0.70196078;
%  0.74509804, 0.72941176, 0.85490196;
%  0.98431373, 0.50196078, 0.44705882;
%  0.50196078, 0.69411765, 0.82745098;
%  0.99215686, 0.70588235, 0.38431373;
%  0.70196078, 0.87058824, 0.41176471;
%  0.98823529, 0.80392157, 0.89803922;
%  0.85098039, 0.85098039, 0.85098039;
%  0.7372549,  0.50196078, 0.74117647;
%  0.8,        0.92156863, 0.77254902;
%  1.,         0.92941176, 0.43529412] % Set3
 
cmap =  [0.12156863, 0.46666667, 0.70588235;
 1.,        0.49803922, 0.05490196;
 0.17254902, 0.62745098, 0.17254902;
 0.83921569, 0.15294118, 0.15686275;
 0.58039216, 0.40392157, 0.74117647;
 0.54901961, 0.3372549,  0.29411765;
 0.89019608, 0.46666667, 0.76078431;
 0.49803922, 0.49803922, 0.49803922;
 0.7372549,  0.74117647, 0.13333333;
 0.09019608, 0.74509804, 0.81176471]; %tab10 

figure('units','normalized','outerposition',[0.2 0.2 0.15 0.35],'Color','w');
for iSub = 1: length(subjs)
 plot(agg_impv(agg_impv(:, 3) == iSub, 1:2)', 'o-', 'LineWidth', 1.5, 'Color', cmap(iSub,:), 'MarkerFaceColor', cmap(iSub,:),'MarkerSize', 8); hold on;     
end 
xlim([0.5 2.5]);  xticks([1 2]); 
ylim([0 50]); yticks([0: 10: 50]); 
xticklabels({'round1', 'round2'}); 
set(gca,'fontsize', 25); 
ylabel('number of moves'); 
% title(num2str(mean(ratio(:,1)))); 

figure('units','normalized','outerposition',[0.2 0.2 0.15 0.35],'Color','w');
for iSub = 1: length(subjs)
 plot(agg_same(agg_same(:, 3) == iSub, 1:2)', 'o-', 'LineWidth', 1.5, 'Color', cmap(iSub,:), 'MarkerFaceColor', cmap(iSub,:),'MarkerSize', 8); hold on;     
end 
plot(agg_same(agg_same(:, 3) == 5, 1:2)', 'o-', 'LineWidth', 1.5, 'Color', cmap(5,:), 'MarkerFaceColor', cmap(5,:),'MarkerSize', 8); hold on;     
plot(agg_same(agg_same(:, 3) == 2, 1:2)', 'o-', 'LineWidth', 1.5, 'Color', cmap(2,:), 'MarkerFaceColor', cmap(2,:),'MarkerSize', 8); hold on;     
plot(agg_same(agg_same(:, 3) == 8, 1:2)', 'o-', 'LineWidth', 1.5, 'Color', cmap(8,:), 'MarkerFaceColor', cmap(8,:),'MarkerSize', 8); hold on;     
plot(agg_same(agg_same(:, 3) == 7, 1:2)', 'o-', 'LineWidth', 1.5, 'Color', cmap(7,:), 'MarkerFaceColor', cmap(7,:),'MarkerSize', 8); hold on;     
plot(agg_same(agg_same(:, 3) == 4, 1:2)', 'o-', 'LineWidth', 1.5, 'Color', cmap(4,:), 'MarkerFaceColor', cmap(4,:),'MarkerSize', 8); hold on;     
% plot(agg_same(agg_same(:, 3) == 6, 1:2)', 'o-', 'LineWidth', 1.5, 'Color', cmap(6,:), 'MarkerFaceColor', cmap(6,:),'MarkerSize', 8); hold on;     
xlim([0.5 2.5]);  xticks([1 2]); 
ylim([0 50]); yticks([0: 10: 50]); 
xticklabels({'round1', 'round2'}); 
set(gca,'fontsize', 25); 
ylabel('number of moves'); 
% title(num2str(mean(ratio(:,2)))); 

figure('units','normalized','outerposition',[0.2 0.2 0.15 0.35],'Color','w');
for iSub = 1: length(subjs)
 plot(agg_worse(agg_worse(:, 3) == iSub, 1:2)', 'o-', 'LineWidth', 1.5, 'Color', cmap(iSub,:), 'MarkerFaceColor', cmap(iSub,:),'MarkerSize', 8); hold on;     
end 
xlim([0.5 2.5]);  xticks([1 2]); 
ylim([0 50]); yticks([0: 10: 50]); 
xticklabels({'round1', 'round2'}); 
set(gca,'fontsize', 25); 
ylabel('number of moves'); 
% title(num2str(mean(ratio(:,3)))); 

% MGS - unified figure for all subjects
ratio_impv(:,1) = sum(ratio(:,1:2)'); % pull together less steps and equal #
ratio_impv(:,2) = ratio(:,3); % more steps
scatter_loc = 1 + 0.2*randn(1,10);
figure('units','normalized','outerposition',[0.2 0.2 0.15 0.35],'Color','w');
hold all
bar(1,100*mean(ratio_impv(:,1)),'facecolor',[0.8 0.8 0.8])
errorbar(1,100*mean(ratio_impv(:,1)),100*std(ratio_impv(:,1))/sqrt(length(subjs)),'.k');
for iSub = 1: length(subjs)
    plot(scatter_loc(iSub), 100*ratio_impv(iSub,1), '.', 'Color', cmap(iSub,:), 'MarkerFaceColor', cmap(iSub,:),'MarkerSize', 20); hold on;     
end 
xlim([0 2]);  xticks(''); 
ylim([0 120]); yticks(100*[0 0.5 1]); 
plot(get(gca,'xlim'),100*[0.5 0.5],'k')
% xticklabels({'impv', 'worst'}); 
set(gca,'fontsize', 12); 
ylabel('Percentage (%)'); 
% title(num2str(mean(ratio(:,3)))); 
savename = ['MGS_ratioImpv_n', num2str(length(subjs))];
saveas(gcf,fullfile(dirs.banal, [savename '.tiff']),'tiff'); 



%% proportion of mazes that people do equally well or better in b2 
addpath(genpath('/Users/kaminkim/Documents/codes/useful'));
Gpro = nan(length(subjs), 1); 
savename = ['ratioNMimpv_n', num2str(length(subjs))];
newA4figure(savename); 
axes('position',[0.1 0.1 0.1 0.2])
for iSub = 1: length(subjs)
    subj = subjs{iSub}; 
    load (fullfile([dirs.banal 'nav_impv_v3/'], [subj '_maze_nav_impv.mat']));    
    Gpro(iSub) = sum(maze_nav_impv.nMoves_impv)/height(maze_nav_impv);  
    plot([1+(0.5-rand(1))/15], [Gpro(iSub)], 'bo', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'w', 'MarkerSize', 4); hold on;
end 
p = violin(Gpro, 'mc', [], 'medc', 'k') ;
p.FaceColor = [0.7 0.7 0.7];
p.FaceAlpha = 0.3;
xlim([0.2 1.8]); ylim([0 1.3]); 
xticks([]); yticks([0 0.5 1]);
set(gca,'TickLength',[0.05, 0.01])
set(gca, 'FontSize', 12); box off
saveas(gcf,fullfile(dirs.banal, [savename '.jpg']),'jpg'); 

%% plot mean RT per event types
addpath(genpath('/Users/kaminkim/Documents/codes/useful'));
for iSub = 1: length(subjs)
    subj = subjs{iSub}; 
    load([dirs.banal,'sq_table/' subj '_sq_table.mat']); 
    
    sq_table(sq_table.block == 3,:) = []; 
    [~, ~, rm_bool] = trim_mean(sq_table.rt , 3);        
    
    if iSub == 1, EvTypes = unique(sq_table.square);    end     
    for iType = 1: length(EvTypes)      
        tmp_rt = sq_table.rt(strcmp(sq_table.square, EvTypes{iType}) & ~rm_bool, :);
        mRT.(EvTypes{iType})(iSub) = mean(tmp_rt);
        clear tpm_rt; 
    end
    
    tmp_ridx = find(strcmp(sq_table.square, 'G'));
    tmp_rt = sq_table.rt(tmp_ridx-1, :);
    mRT.preG(iSub) = trim_mean(tmp_rt, 3); clear tmp_ridx tmp_rt;   
       
    tmp_ridx = find(strcmp(sq_table.square, 'O'));
    tmp_rt = sq_table.rt(tmp_ridx-1, :);
    mRT.preO(iSub) = trim_mean(tmp_rt, 3); clear tmp_ridx tmp_rt;   

    tmp_ridx = find(strcmp(sq_table.square, 'D'));
    tmp_rt = sq_table.rt(tmp_ridx-1, :);
    mRT.preD(iSub) = trim_mean(tmp_rt, 3); clear tmp_ridx tmp_rt;   
end
% 
% % adding blocks 
% for iSub = 1: length(subjs)
%     subj = subjs{iSub}; 
%     load([dirs.banal,'sq_table/' subj '_sq_table.mat']); 
%     
%     if iSub == 1, EvTypes = unique(sq_table.square);    end     
% 
%     ridx_b1 = find(sq_table.block == 1);
%     ridx_b2 = find(sq_table.block == 2);    
%     
%     tmp_ridx = intersect(find(strcmp(sq_table.square, 'G')), ridx_b1);
%     mRT.preG_b1(iSub) = mean(sq_table.rt(tmp_ridx-1, :));       
%     tmp_ridx = intersect(find(strcmp(sq_table.square, 'O')), ridx_b1);
%     mRT.preO_b1(iSub) = mean(sq_table.rt(tmp_ridx-1, :));       
%     tmp_ridx = intersect(find(strcmp(sq_table.square, 'D')), ridx_b1);
%     mRT.preD_b1(iSub) = mean(sq_table.rt(tmp_ridx-1, :));       
% 
%     tmp_ridx = intersect(find(strcmp(sq_table.square, 'G')), ridx_b2);
%     mRT.preG_b2(iSub) = mean(sq_table.rt(tmp_ridx-1, :));       
%     tmp_ridx = intersect(find(strcmp(sq_table.square, 'O')), ridx_b2);
%     mRT.preO_b2(iSub) = mean(sq_table.rt(tmp_ridx-1, :));       
%     tmp_ridx = intersect(find(strcmp(sq_table.square, 'D')), ridx_b2);
%     mRT.preD_b2(iSub) = mean(sq_table.rt(tmp_ridx-1, :));       
%     
%     tmp_ridx = intersect(find(strcmp(sq_table.square, 'O')), ridx_b1);
%     mRT.O_b1(iSub) = mean(sq_table.rt(tmp_ridx, :));       
%     tmp_ridx = intersect(find(strcmp(sq_table.square, 'D')), ridx_b1);
%     mRT.D_b1(iSub) = mean(sq_table.rt(tmp_ridx, :));       
%     tmp_ridx = intersect(find(strcmp(sq_table.square, 'G')), ridx_b1);
%     mRT.G_b1(iSub) = mean(sq_table.rt(tmp_ridx, :));       
%     tmp_ridx = intersect(find(strcmp(sq_table.square, 'X')), ridx_b1);
%     mRT.X_b1(iSub) = mean(sq_table.rt(tmp_ridx, :));    
%     
%     tmp_ridx = intersect(find(strcmp(sq_table.square, 'O')), ridx_b2);
%     mRT.O_b2(iSub) = mean(sq_table.rt(tmp_ridx, :));       
%     tmp_ridx = intersect(find(strcmp(sq_table.square, 'D')), ridx_b2);
%     mRT.D_b2(iSub) = mean(sq_table.rt(tmp_ridx, :));       
%     tmp_ridx = intersect(find(strcmp(sq_table.square, 'G')), ridx_b2);
%     mRT.G_b2(iSub) = mean(sq_table.rt(tmp_ridx, :));       
%     tmp_ridx = intersect(find(strcmp(sq_table.square, 'X')), ridx_b2);
%     mRT.X_b2(iSub) = mean(sq_table.rt(tmp_ridx, :));    
%     clear ridx* tmp_ridx 
% end

% New organization 
% rtcell = {mRT.X mRT.N mRT.D mRT.E mRT.O mRT.G, mRT.preG, mRT.preO};
% mnames = {'new maze', 'open path', 'choice point', 'goal entering', 'path opened', 'goal revealed', 'expect\npath opening', 'expect\ngoal revealing'}; 
% figure; 
% p = violin(rtcell, 'mc', [], 'medc', 'k') ;
% for ii = 1: length(mnames)
%     for jj = 1: length(rtcell{ii}), shake_by(jj, 1) = rand(1); end 
%     p(ii).FaceColor = [0.7 0.7 0.7]; p(ii).FaceAlpha = 0.3; p(ii).EdgeColor = [1 1 1];  
%     plot([ii+(0.5-shake_by)/8], rtcell{ii}, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 7); hold on; legend off; 
%     clear shake_by 
% end 
% % xlim([0 10]); ylim([-500 5000]); 
% xticks([1:9]); xticklabels(mnames);  
% % yticks([-500: 500: 5000]);
% set(gca, 'FontSize', 16); box off
% ylabel('RT(ms)'); xlabel('move type'); 
% saveas(gcf,fullfile(dirs.banal, ['RT_perEvType_n' num2str(length(subjs)) '_trimmed_mean_organized_forpaper_wPre_b12.tiff']),'tiff'); 

% New organization2: for MS  
rtcell = {mRT.X, mRT.D, mRT.G, mRT.E, mRT.N };
mnames = {'new maze',  'choice point', 'goal revealed', 'goal entering', 'open path'}; 
figure; 
p = violin(rtcell, 'mc', [], 'medc', 'k') ;
for ii = 1: length(mnames)
    for jj = 1: length(rtcell{ii}), shake_by(jj, 1) = rand(1); end 
    p(ii).FaceColor = [0.7 0.7 0.7]; p(ii).FaceAlpha = 0.3; p(ii).EdgeColor = [1 1 1];  
    plot([ii+(0.5-shake_by)/8], rtcell{ii}, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5); hold on; legend off; 
    clear shake_by 
end 
% xlim([0 10]); ylim([-500 5000]); 
% xticks([1:9]); xticklabels(mnames);  
% yticks([-500: 500: 5000]);
set(gca, 'FontSize', 25); box off
ylabel('RT(ms)'); xlabel('move type'); 
saveas(gcf,fullfile(dirs.banal, ['RT_perEvType_n' num2str(length(subjs)) '_trimmed_mean_organized_forpaper_woPre_b12.tiff']),'tiff'); 

%% New organization3: splitting B1 and B2 for MS 
cmapGr = brewermap(6,'Greys');
cmap(1,:) = cmapGr(3,:);
cmap(2,:) = cmapGr(5,:);

clear mRT
for iSub = 1: length(subjs)
    subj = subjs{iSub}; 
    load([dirs.banal,'sq_table/' subj '_sq_table.mat']); 
    
    sq_table(sq_table.block == 3,:) = []; 
    [~, ~, rm_bool] = trim_mean(sq_table.rt , 3);        
    
    if iSub == 1, EvTypes = unique(sq_table.square);    end
    for iType = 1: length(EvTypes)
        for iB = 1:2
            tmp_rt = sq_table.rt(strcmp(sq_table.square, EvTypes{iType}) & (sq_table.block == iB)& (sq_table.opt_path_impv12 == 1) & ~rm_bool, :);
            mRT.(EvTypes{iType})(iSub, iB) = mean(tmp_rt);
            clear tpm_rt;
        end
    end
    
    for iB = 1:2
        tmp_ridx = find(strcmp(sq_table.square, 'G') & (sq_table.block == iB) & (sq_table.opt_path_impv12 == 1) & ~rm_bool);
        tmp_rt = sq_table.rt(tmp_ridx-1, :);
        mRT.preG(iSub, iB) = trim_mean(tmp_rt, 3); clear tmp_ridx tmp_rt;

        tmp_ridx = find(strcmp(sq_table.square, 'O') & (sq_table.block == iB));
        tmp_rt = sq_table.rt(tmp_ridx-1, :);
        mRT.preO(iSub, iB) = trim_mean(tmp_rt, 3); clear tmp_ridx tmp_rt;

        tmp_ridx = find(strcmp(sq_table.square, 'D') & (sq_table.block == iB));
        tmp_rt = sq_table.rt(tmp_ridx-1, :);
        mRT.preD(iSub, iB) = trim_mean(tmp_rt, 3); clear tmp_ridx tmp_rt;
    end
end
for iType = 1:length(EvTypes)
    [pV(iType),h] = signrank(mRT.(EvTypes{iType})(:,1),mRT.(EvTypes{iType})(:,2));
end

newA4figure('')
rtcell = {mRT.X(:,1), mRT.X(:,2)};
mnames = {'Outset'}; 
axes('position',[0.1 0.1 0.1 0.4])
p = violin(rtcell, 'mc', [], 'medc', 'k') ;
for ii = 1: length(rtcell)
    for jj = 1: length(rtcell{ii}), shake_by(jj, 1) = rand(1); end 
    if mod(ii,2)
        p(ii).FaceColor = cmap(1,:); p(ii).FaceAlpha = 0.3; p(ii).EdgeColor = [1 1 1];  
    else
        p(ii).FaceColor = cmap(2,:); p(ii).FaceAlpha = 0.3; p(ii).EdgeColor = [1 1 1];  
    end
    plot([ii+(0.5-shake_by)/8], rtcell{ii}, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5); hold on; legend off; 
    clear shake_by 
end 
set(gca,'Xtick',[1.5 ],  'XTickLabel',mnames{1},'XTickLabelRotation',0)     
ylabel('RT(ms)'); 
set(gca, 'ylim',[0 3000],'ytick',[0 1500 3000])

axes('position',[0.28 0.1 0.2 0.4])
rtcell = {mRT.preG(:,1), mRT.preG(:,2), mRT.N(:,1), mRT.N(:,2)};
mnames = { 'Pre-Goal','Open path'}; 
p = violin(rtcell, 'mc', [], 'medc', 'k') ;
for ii = 1: length(rtcell)
    for jj = 1: length(rtcell{ii}), shake_by(jj, 1) = rand(1); end 
    if mod(ii,2)
        p(ii).FaceColor = [0.8 0.8 0.8]; p(ii).FaceAlpha = 0.3; p(ii).EdgeColor = [1 1 1];  
    else
        p(ii).FaceColor = [0.4 0.4 0.4]; p(ii).FaceAlpha = 0.3; p(ii).EdgeColor = [1 1 1];  
    end
    plot([ii+(0.5-shake_by)/8], rtcell{ii}, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5); hold on; legend off; 
    clear shake_by 
end 
set(gca,'Xtick',[1.5 3.5],  'XTickLabel',mnames,'XTickLabelRotation',0) 
set(gca, 'ylim',[0 800],'ytick',[0 400 800])

saveas(gcf,fullfile(dirs.banal, ['RT_perEvType_n' num2str(length(subjs)) '_trimmed_mean_organized_forpaper_O_G_N_impv_b12.tiff']),'tiff'); 

%%
% New organization3: splitting B1 and B2 for MS 
for ii_b = 1:2 % compare improved and degraded trials across blocks

    clear mRT
    for iSub = 1: length(subjs)
        subj = subjs{iSub};
        load([dirs.banal,'sq_table/' subj '_sq_table.mat']);

        sq_table(sq_table.block == 3,:) = [];
        [~, ~, rm_bool] = trim_mean(sq_table.rt , 3);

        if iSub == 1, EvTypes = unique(sq_table.square);    end
        for iType = 1: length(EvTypes)
            for iB = 1:2
                if iB == 1
                    tmp_rt = sq_table.rt(strcmp(sq_table.square, EvTypes{iType}) & (sq_table.block == ii_b)& (sq_table.opt_path_impv12 == 1) & ~rm_bool, :);
                else
                    tmp_rt = sq_table.rt(strcmp(sq_table.square, EvTypes{iType}) & (sq_table.block == ii_b)& (sq_table.opt_path_impv12 == 0) & ~rm_bool, :);
                end
                mRT.(EvTypes{iType})(iSub, iB) = mean(tmp_rt);
                clear tpm_rt;
            end
        end

        for iB = 1:2
            if iB == 1
                tmp_ridx = find(strcmp(sq_table.square, 'G') & (sq_table.block == 2) & (sq_table.opt_path_impv12 == 1));
            else
                tmp_ridx = find(strcmp(sq_table.square, 'G') & (sq_table.block == 2) & (sq_table.opt_path_impv12 == 0));
            end
            tmp_rt = sq_table.rt(tmp_ridx-1, :);
            mRT.preG(iSub, iB) = trim_mean(tmp_rt, 3); clear tmp_ridx tmp_rt;

            tmp_ridx = find(strcmp(sq_table.square, 'O') & (sq_table.block == ii_b));
            tmp_rt = sq_table.rt(tmp_ridx-1, :);
            mRT.preO(iSub, iB) = trim_mean(tmp_rt, 3); clear tmp_ridx tmp_rt;

            tmp_ridx = find(strcmp(sq_table.square, 'D') & (sq_table.block == ii_b));
            tmp_rt = sq_table.rt(tmp_ridx-1, :);
            mRT.preD(iSub, iB) = trim_mean(tmp_rt, 3); clear tmp_ridx tmp_rt;
        end
    end

    for iType = 1:length(EvTypes)
        [pV(iType,ii_b),h] = signrank(mRT.(EvTypes{iType})(:,1),mRT.(EvTypes{iType})(:,2));
    end

    figname = sprintf('RT_perEvType_n%d_trimmed_mean_organized_forpaper_O_G_N_impvDeg_b%d',length(subjs),ii_b);
    newA4figure(figname)
    cmap = brewermap(2,'Dark2');
    rtcell = {mRT.X(:,1), mRT.X(:,2)};
    mnames = {'Outset'};
    axes('position',[0.1 0.1 0.1 0.4])
    p = violin(rtcell, 'mc', [], 'medc', 'k') ;
    for ii = 1: length(rtcell)
        for jj = 1: length(rtcell{ii}), shake_by(jj, 1) = rand(1); end
        if mod(ii,2)
            p(ii).FaceColor = cmap(1,:); p(ii).FaceAlpha = 0.3; p(ii).EdgeColor = cmap(1,:);
        else
            p(ii).FaceColor = cmap(2,:); p(ii).FaceAlpha = 0.3; p(ii).EdgeColor = cmap(2,:);
        end
        plot([ii+(0.5-shake_by)/8], rtcell{ii}, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5); hold on; legend off;
        clear shake_by
    end
    set(gca,'Xtick',[1.5 ],  'XTickLabel',mnames{1},'XTickLabelRotation',0)
    ylabel('RT(ms)');
    set(gca, 'ylim',[0 3500],'ytick',[0 1500 3500])

    axes('position',[0.28 0.1 0.2 0.4])
    rtcell = {mRT.preG(:,1), mRT.preG(:,2), mRT.N(:,1), mRT.N(:,2)};
    mnames = { 'Pre-Goal','Open path'};
    p = violin(rtcell, 'mc', [], 'medc', 'k') ;
    for ii = 1: length(rtcell)
        for jj = 1: length(rtcell{ii}), shake_by(jj, 1) = rand(1); end
        if mod(ii,2)
            p(ii).FaceColor = cmap(1,:); p(ii).FaceAlpha = 0.3; p(ii).EdgeColor = cmap(1,:);
        else
            p(ii).FaceColor = cmap(2,:); p(ii).FaceAlpha = 0.3; p(ii).EdgeColor = cmap(2,:);
        end
        plot([ii+(0.5-shake_by)/8], rtcell{ii}, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5); hold on; legend off;
        clear shake_by
    end
    set(gca,'Xtick',[1.5 3.5],  'XTickLabel',mnames,'XTickLabelRotation',0)
    set(gca, 'ylim',[0 800],'ytick',[0 400 800])

    saveas(gcf,fullfile(dirs.banal, figname),'tiff');

end


%% mixed-model stats on RT: trimmed mean & mean-centered
% model: RT ~ ttype * blk + (1|subj) 
 
% initialize variable structs - mega, all trial types 
data_all_ttype = {'X', 'N', 'D', 'E', 'O', 'G', 'preO', 'preG'}; 
clear agg_data
agg_data.RT = [];
agg_data.ttype = [];
agg_data.blk = [];
agg_data.subj = [];
agg_data.impv = [];
% fill in variable structs 
for iSub = 1: length(subjs)
    subj = subjs{iSub}; 
    load([dirs.banal,'sq_table/' subj '_sq_table.mat']); 
    [~, ~, rm_bool] = trim_mean(sq_table.rt , 3);  
    sq_table.rt =  (sq_table.rt-nanmean(sq_table.rt))/nanstd(sq_table.rt);
    
    for iType = 1: length(data_all_ttype)
        ridx = contains(sq_table.square, data_all_ttype{iType});%  & sq_table.opt_path_impv12 == 1;
        ridx = find(ridx & ~rm_bool);  
        agg_data.RT = cat(1, agg_data.RT, sq_table.rt(ridx));
        agg_data.ttype = cat(1, agg_data.ttype, repmat([iType], [length(ridx), 1]));
        agg_data.blk =  cat(1, agg_data.blk, sq_table.block(ridx));
        agg_data.impv =  cat(1, agg_data.impv, sq_table.opt_path_impv12(ridx));
        agg_data.subj =  cat(1, agg_data.subj, repmat([iSub],  [length(ridx), 1]));   

        clear ridx 
    end
end
% look at summary stats
rt_table = table(agg_data.ttype, agg_data.impv, agg_data.RT, agg_data.blk, agg_data.subj, 'VariableNames', {'type', 'impv','rt', 'block', 'subj'});
rt_table(rt_table.block == 3, :) = [];
grptab = grpstats(rt_table, {'type', 'block'}, {'mean', 'std'});
model = fitlme(rt_table, 'rt ~ impv + impv*block + (1|subj)');

rt_table_impv = rt_table(rt_table.impv == 1,:);
model = fitlme(rt_table_impv, 'rt ~ type + block + type*block + (1|subj)');

rt_table_deg = rt_table(rt_table.impv == 0,:);
model = fitlme(rt_table_deg, 'rt ~ type + block + type*block + (1|subj)');

rt_table_b1 = rt_table(rt_table.block == 1,:);
model = fitlme(rt_table_b1, 'rt ~ type + impv + type*impv + (1|subj)');
rt_table_b2 = rt_table(rt_table.block == 2,:);
model = fitlme(rt_table_b2, 'rt ~ type + impv + type*impv + (1|subj)');

model
m = [];
m(1, 1) = grptab.mean_rt(grptab.type == 3 & grptab.block == 1);
m(1, 2) = grptab.mean_rt(grptab.type == 3 & grptab.block == 2);
% m(1, 3) = grptab.mean_rt(grptab.type == 2 & grptab.block == 3);
m(2, 1) = grptab.mean_rt(grptab.type == 2 & grptab.block == 1);
m(2, 2) = grptab.mean_rt(grptab.type == 2 & grptab.block == 2);
%m(2, 3) = grptab.mean_rt(grptab.type == 3 & grptab.block == 3);

se = [];
se (1, 1) = grptab.std_rt(grptab.type == 3 & grptab.block == 1)/sqrt(1);
se (1, 2) = grptab.std_rt(grptab.type == 3 & grptab.block == 2)/sqrt(9);
%se (1, 3) = grptab.std_rt(grptab.type == 2 & grptab.block == 3)/sqrt(9);
se (2, 1) = grptab.std_rt(grptab.type == 2 & grptab.block == 1)/sqrt(9);
se (2, 2) = grptab.std_rt(grptab.type == 2 & grptab.block == 2)/sqrt(9);
%se (2, 3) = grptab.std_rt(grptab.type == 3 & grptab.block == 3)/sqrt(9);

figure; 
hold all
errorbar([0.95 1.95], m(1, :)', se(1, :)',  'Color', [0.4 0.4 0.4], 'LineWidth', 1.5); hold on; 
errorbar([1.05 2.05], m(2, :)', se(2, :)',  'Color',  cmap(1, :), 'LineWidth', 1.5); 
xlim([0.5 2.5]); ylim([-1 1.5]); 
xticks([1 2 ]); yticks([-1:0.5:1]);
box off; 
xlabel('block'); ylabel('z (RT)')
set(gca, 'FontSize', 16); 
legend({'Choice points', 'Open Path'});
legend boxoff; 
%%

% set up params and data for stats
% nperm = 500; 
contrasts = {'N_X', 'N_D', 'N_E', 'N_O', 'N_G', 'O_G', 'N_preG'}; 
pmat = nan(4, length(contrasts)); 
for iCon = 1: length(contrasts) 
    tmpstr = contrasts{iCon}; sep = strfind(tmpstr, '_');
    conds = [find(strcmp(data_all_ttype, tmpstr(1: sep-1))), find(strcmp(data_all_ttype, tmpstr(sep+1: end)))]; 
    Y = cat(1, agg_data.RT(agg_data.ttype == conds(1)), agg_data.RT(agg_data.ttype == conds(2))); 
    X1 = cat(1, agg_data.ttype(agg_data.ttype == conds(1)), agg_data.ttype(agg_data.ttype == conds(2))); 
    X2 = cat(1, agg_data.blk(agg_data.ttype == conds(1)), agg_data.blk(agg_data.ttype == conds(2))); 
    X3 = cat(1, agg_data.subj(agg_data.ttype == conds(1)), agg_data.subj(agg_data.ttype == conds(2))); 
    Data = table(Y, X1, X2, X3, 'VariableNames', {'RT', 'ttype', 'block', 'subj'}); 
    
    % if limiting to blks 1&2
%     Data(Data.block == 3, :) = [];
    
    model1.(contrasts{iCon}) = fitlme(Data, 'RT ~ ttype*block + (1|subj)');
    [beta, betanames, stats]  = fixedEffects(model1.(contrasts{iCon})); 
    pmat(:, iCon) = stats.pValue;    
end 

% save p-vals
clear p;
p{1,1} = 'intercept'; 
p{2,1} = 'ttype'; 
p{3,1} = 'block'; 
p{4,1} = 'intx'; 
for iCon = 1: length(contrasts) 
    for ii = 1: 4
        p{ii, iCon+1} = pmat(ii, iCon);  
    end
end 
header = cat(2, {'vars'}, contrasts); 
p = cat(1, header, p); 
cell2csv(fullfile(dirs.banal, ['zRT_mm_contrast_impv_p_n', num2str(length(subjs)), '_wPre.csv']), p);   



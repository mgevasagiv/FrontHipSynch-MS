% MAZE12_OFC_Gamma_stats
% Factorial stats analyzing TFR diagrams performed by MAZE10_calc_TFR_event

% Separate TFR plots per OFC area
clear all; close all

%% Define lateral and medial contacts
% manuscript - 19, 20, 23
lOFC_contacts = [11  31  32  61];

% manuscript - 21, 24, 25
mOFC_contacts = [41  42  71  81];

cmap = brewermap(3,'PRGn');
cmap(2,:) = 0.8*cmap(2,:);
for ii = 1:3
    cmap_behav{ii} = cmap(ii,:);
end
% 'X' , 'N', 'G'
cmap_behav{4} = cmap_behav{3}*0.6; % 'preG'
%%
% Reorganize colors for the needed order - 'x' 'g' 'preg' 'n'
cmap_behav_g{1} = cmap_behav{1};
cmap_behav_g{2} = cmap_behav{3};
cmap_behav_g{3} = cmap_behav{4};
cmap_behav_g{4} = cmap_behav{2};


%%
% Collecting and Plotting the output of
addpath(genpath('E:\Dropbox\Code\Useful_code'))
addpath(genpath('C:\Users\mgeva\Documents\GitHub\miscTools\distributionPlot'))

srate = 512;
c_count = 1;
agg_data = {};
fb_theta = [4 8];
fb_alpha = [9 11];
fb_alphaBeta = [9 20];
fb_lgamma = [30 50];
fb_mgamma = [50 100] ;
win(1,:) = [0.2 0.35];
win(2,:) = [-0.25 -0.15];
win(3,:) = [0.15 0.25];
minEvents = 3;

roi_list = {'HPC','OFC'};
% Full list - subjects = {'s01', 's06', 's08', 's09', 's12', 's14', 's15', 's16', 's17', 's18'};
evlist = [{'X'}    {'D'}    {'G'}    {'E'}    {'O'}    {'N'}    {'preG'}    {'preO'}    {'null'}];

savedirs = {'E:\MAZE_1.0\analysis\TF\TFR_grp\HPC\population',...
            'E:\MAZE_1.0\analysis\TF\TFR_grp\OFC\population'};

plotdir = 'E:\MAZE_1.0\analysis\TF\TFR_grp\population\';

%%
cond = [0 1 3 4];
text_cond = {'impvTrial_B1','impvTrial_B2','detTrial_B1','detTrial_B2'};

%%
roi_i = 2; cond = 0;
filename = sprintf('CenteredAvgEvents_allElectrodesSubj_%s_%d.mat',roi_list{roi_i},cond);
disp(fullfile(savedirs{roi_i},filename))
mm = matfile(fullfile(savedirs{roi_i},filename));
subjects = mm.subjects;
evlist = mm.evlist;
params = mm.params;


%% Create a table for start/end behavioral timepoints 
for ii_a = 1:2

    roi_i = 2;
    
    if ii_a == 1
        cond = 1; % BLK 2 for impv routes
        blk = 2;

    else
        cond = 0; % BLK 1 for impv routes
        blk = 1;
    end
    
    filename = sprintf('CenteredAvgEvents_allElectrodesSubj_%s_%d.mat',roi_list{roi_i},cond);
    disp(fullfile(savedirs{roi_i},filename))
    mm = matfile(fullfile(savedirs{roi_i},filename));
        
    % GammaBand_trial = average over TFR tXf matrix  
    GammaBehavTable = mm.GammaBehavTable;
    A = GammaBehavTable.('N');
    B = GammaBehavTable.('X');
    C = GammaBehavTable.('G');
    D = GammaBehavTable.('preG');

    clear label_e
    for i = 1:height(A)
        label_e{i} = 'N';
    end
    A.event = label_e';
    clear label_e

    for i = 1:height(B)
        label_e{i} = 'X';
    end
    B.event = label_e';
    clear label_e
    for i = 1:height(C)
        label_e{i} = 'G';
    end
    C.event = label_e';
    clear label_e
    for i = 1:height(D)
        label_e{i} = 'preG';
    end
    D.event = label_e';

    T = [A; B; C; D];

    T.blk = ones(height(T),1)*blk;

    if ii_a == 1
        T_1 = T;
        clear T
    end
    if ii_a == 2
        T = [T; T_1];
    end

end

% statistics when aggregating all start/end timepoints
Data = T;
formula = 'Gamma ~ event + blk + event*blk + (1|subj)';
lme1 = fitlme(Data,formula);
anova(lme1)
betaTable = lme1.Coefficients;

events = {'X','G','preG'};
clear ANO
for ii_e = 1:length(events)
    T_sub = T( strcmp(events{ii_e},T.event) |  strcmp('N',T.event),: ) ;
    formula = 'Gamma ~ event + blk + event*blk + (1|subj)';
    lme1 = fitlme(T_sub,formula);
    ANO{ii_e} = anova(lme1);
end

%%
bootstrapN = 1000;

fig_name = 'GammaPower_behav_block';
newA4figure('GammaPower_behav_block')
for ii_e = 1:length(events)
    clear stringT
    data_x = T.Gamma( strcmp(events{ii_e},T.event) ,: ) ;
    data_y = T.Gamma( strcmp('N',T.event) ,: ) ;

    effect = meanEffectSize(data_x,data_y, Effect ="cohen",ConfidenceIntervalType="bootstrap", ...
        BootstrapOptions=statset(UseParallel=true),NumBootstraps= bootstrapN);
    hold all
    subplot(4,3,ii_e)
    h = gardnerAltmanPlot(data_x,data_y,Paired=false,Effect="robustcohen",Alpha=0.05);
    stringT{1} = sprintf('Main effect (behav)');
    stringT{2} = sprintf('P = %2.2e, Robust Cohens d = %2.2d',ANO{ii_e}.pValue(2),effect.Effect);
    title(stringT)
    if ii_e == 1
        ylabel('OFC Gamma power');
    end
    ax = get(gca);
    set(gca,'xticklabel',{events{ii_e},'N', 'Robust Cohen''s d'});
    h(1).CData = cmap_behav_g{ii_e};
    h(2).CData = cmap_behav_g{4};
end

for ii_e = 1:length(events)
    clear stringT
    data_x = T.Gamma( T.blk ==1 & strcmp(events{ii_e},T.event) ,: ) ;
    data_y = T.Gamma( T.blk ==2 & strcmp(events{ii_e},T.event) ,: ) ;

    effect = meanEffectSize(data_x,data_y, Effect ="cohen",ConfidenceIntervalType="bootstrap", ...
        BootstrapOptions=statset(UseParallel=true),NumBootstraps=bootstrapN);
    hold all
    subplot(4,3,3 + ii_e)
    h = gardnerAltmanPlot(data_x,data_y,Paired=false,Effect="robustcohen",Alpha=0.05);
    stringT{1} = sprintf('Main effect (blk)');
    stringT{1} = sprintf('P = %2.2e, Robust Cohens d = %2.2d',ANO{ii_e}.pValue(3), effect.Effect);
    title(stringT)
    if ii_e == 1
        ylabel('OFC Gamma power');
    end
    ax = get(gca);
    set(gca,'xticklabel',{'Blk 1','Blk 2', 'Robust Cohen''s d'});
     h(1).CData = cmap_behav_g{ii_e};
     h(2).CData = cmap_behav_g{4};
end

for ii_e = 1:length(events)
    clear stringT
    data_x = T.Gamma( T.blk ==1 & strcmp(events{ii_e},T.event) ,: ) ;
    data_y = T.Gamma( T.blk ==1 & strcmp('N',T.event) ,: ) ;

    effect = meanEffectSize(data_x,data_y, Effect ="cohen",ConfidenceIntervalType="bootstrap", ...
        BootstrapOptions=statset(UseParallel=true),NumBootstraps=bootstrapN);
    hold all
    subplot(4,3,6 + ii_e)
    h = gardnerAltmanPlot(data_x,data_y,Paired=false,Effect="robustcohen",Alpha=0.05);
    stringT = sprintf('Blk1, Robust Cohens d = %2.2d',effect.Effect);
    title(stringT)
    if ii_e == 1
        ylabel('OFC Gamma power');
    end    
         h(1).CData = cmap_behav_g{ii_e};
     h(2).CData = cmap_behav_g{4};

    ax = get(gca);
    set(gca,'xticklabel',{events{ii_e},'N',  'Robust Cohen''s d'})
    XLIM = get(gca,'xlim');
    YLIM = get(gca,'ylim');

end

for ii_e = 1:length(events)
    clear stringT
    data_x = T.Gamma( T.blk == 2 & strcmp(events{ii_e},T.event) ,: ) ;
    data_y = T.Gamma( T.blk == 2 & strcmp('N',T.event) ,: ) ;

    effect = meanEffectSize(data_x,data_y, Effect ="cohen",ConfidenceIntervalType="bootstrap", ...
        BootstrapOptions=statset(UseParallel=true),NumBootstraps=bootstrapN);
    hold all
    subplot(4,3,9 + ii_e)
    h = gardnerAltmanPlot(data_x,data_y,Paired=false,Effect="robustcohen",Alpha=0.05);
    stringT = sprintf('Blk2, Robust Cohens d = %2.2d', effect.Effect);
    title(stringT)
         h(1).CData = cmap_behav_g{ii_e};
     h(2).CData = cmap_behav_g{4};
    if ii_e == 1
        ylabel('OFC Gamma power');
    end
    ax = get(gca);
    set(gca,'xticklabel',{events{ii_e},'N',  'Robust Cohen''s d'});
        XLIM = get(gca,'xlim');
    YLIM = get(gca,'ylim');

    stringT2{1} = sprintf('Interaction ');
    stringT2{2} = sprintf('blk x location');
    stringT2{3} = sprintf('P = %2.2e',ANO{ii_e}.pValue(4));
    text(XLIM(1)+diff(XLIM)/8,YLIM(1)-diff(YLIM)/3,stringT2);

end
% saveas(gcf, fullfile(plotdir,fig_name),'tiff')


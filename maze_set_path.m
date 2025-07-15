% adding this line just as a test
[ret, name] = system('hostname'); 

if contains(name, 'cns-cdkngw2')
% if contains(name, 'Kamins-iMac.local')
    pj_home_dir = 'E:\MAZE_1.0\' ;
else 
    error ('Your computer is not recognized: check maze_set_path!'); 
end 

dirs.nldata = [pj_home_dir 'data\iEEG\raw_orig\neuralynx\']; 
dirs.rawdata = [pj_home_dir  'data\iEEG\raw\']; 
dirs.prepfile = [pj_home_dir  'data\iEEG\preproc\']; 
dirs.ripple = [pj_home_dir   'data\iEEG\ripple\']; 
dirs.fourier = [pj_home_dir  'data\iEEG\fourier\']; 
dirs.behavioral = [pj_home_dir  'data\behav\']; 
dirs.ctmr = [pj_home_dir  'data\CT_MR\']; 
dirs.elec = [pj_home_dir   'data\electrode_placement\']; 
dirs.ripplerate_ER = [pj_home_dir   'analysis\ripplerate_ER\']; 
dirs.ripplerate_stype = [pj_home_dir   'analysis\ripplerate_ST\']; 
dirs.banal = [pj_home_dir   'analysis\behav\']; 
dirs.TF = [pj_home_dir   'analysis\TF\']; 
dirs.analysis_root = [pj_home_dir   'analysis\']; 
dirs.grp_roi_mask = [pj_home_dir   'analysis\group_emask\']; 
dirnames = fields(dirs); 
for i = 1: length(dirnames)
    if ~exist(dirs.(dirnames{i}))
        mkdir(dirs.(dirnames{i}));
    end 
end 

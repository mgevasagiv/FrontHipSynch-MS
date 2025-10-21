% This function compares navigation performance in block 1 and block 2 - for each maze
% REQUIRES: maze_difficulty.mat and subject's logtable 
% OUTPUT (saves):
% 'maze' table with names (world), associated trial index, maze difficulty, & navigation improvement measures
% 'dpmatch' structs that marks which decision points in block1 and block2
% trials match (1 = match, 0 = no match). 
% v3: matching criteria - D coords 

% ------------------- 2021 June ---------------------------
% ------------ by Radhika Dhanak, Kamin Kim -------

function func_nav_impv_v3 (subj, logdir, outdir) 
    % load maze difficulty and subject's logtable 
    load('maze_difficulty.mat'); 
    load(fullfile(logdir, [subj '_logtable.mat'])); 

    worlds = unique(logtable.world); %finds all unique maze names
    data = nan(30,29);

    %% For each world
    for iW = 1: length(worlds)
        %% Which trial was this maze in? (trialidx within each block)
        mazeTidx_b1 = logtable.tidx(find(strcmp(logtable.world, worlds(iW)) & (logtable.block == 1)));
        mazeTidx_b1 = mazeTidx_b1 (1);
        mazeTidx_b2 = logtable.tidx(find(strcmp(logtable.world, worlds(iW)) & (logtable.block == 2)));
        mazeTidx_b2 = mazeTidx_b2 (1);
        mazeTidx_b3 = logtable.tidx(find(strcmp(logtable.world, worlds(iW)) & (logtable.block == 3)));
        if isempty(mazeTidx_b3) 
            mazeTidx_b3 = 0;
        else
            mazeTidx_b3 = mazeTidx_b3 (1);
        end
        %% Total navigation time for each maze        
        % block 1
        mazeTimestamps_b1 = logtable.timestamp(find(strcmp(logtable.world, worlds(iW)) & (logtable.block == 1)));
        data(iW,2) = mazeTimestamps_b1(end) - mazeTimestamps_b1(1);
        % block 2
        mazeTimestamps_b2 = logtable.timestamp(find(strcmp(logtable.world, worlds(iW)) & (logtable.block == 2)));
        data(iW,3) = mazeTimestamps_b2(end) - mazeTimestamps_b2(1);
        % block 3
        mazeTimestamps_b3 = logtable.timestamp(find(strcmp(logtable.world, worlds(iW)) & (logtable.block == 3)));
        if mazeTidx_b3 == 0 
            data(iW,4) = 0;
        else
            data(iW,4) = mazeTimestamps_b3(end) - mazeTimestamps_b3(1);
        end
        
        % was there any improvement between blocks?    
        if mazeTidx_b3 == 0 
          if data(iW,2) < data(iW,3)
              data(iW,5) = 0;
          else
              data(iW,5) = 1;
          end
        else
            if  data(iW,2) < data(iW,3) || data(iW,2) < data(iW,4)
                data(iW,5) = 0;
            else
                data(iW,5) = 1;
            end
        end
                
        %% Number of moves to complete each maze
        % block 1
        mazeSquaretype_b1 = logtable.squaretype(find(strcmp(logtable.world, worlds(iW)) & (logtable.block == 1)));
        data(iW,6) = length(mazeSquaretype_b1)-1;
        % block 2
        mazeSquaretype_b2 = logtable.squaretype(find(strcmp(logtable.world, worlds(iW)) & (logtable.block == 2)));
        data(iW,7) = length(mazeSquaretype_b2)-1;
        % block 3
        mazeSquaretype_b3 = logtable.squaretype(find(strcmp(logtable.world, worlds(iW)) & (logtable.block == 3)));
        if mazeTidx_b3 == 0
            data(iW,8) = 0;
        else
            data(iW,8) = length(mazeSquaretype_b3)-1;
        end

        % was there any improvement between blocks?    
        if mazeTidx_b3 == 0 
          if data(iW,6)  < data(iW,7)
              data(iW,9) = 0;
          else
              data(iW,9) = 1;
          end
        else
            if  data(iW,6) < data(iW,7) || data(iW,6) < data(iW,8)
                data(iW,9) = 0;
            else
                data(iW,9) = 1;
            end
        end        
        %% Average decision time taken for each maze - limited to decision points matchin in b1&b2
        % For each decition point-pair (b1-b2), compare 
        % blackremains, numsquresopen, and position (maze coordinate)
        dRidx_b1 = find(strcmp(logtable.world, worlds(iW)) & strcmp(logtable.squaretype, 'D') & (logtable.block == 1));
        dRidx_b2 = find(strcmp(logtable.world, worlds(iW)) & strcmp(logtable.squaretype, 'D') & (logtable.block == 2));
        dRidx_b3 = find(strcmp(logtable.world, worlds(iW)) & strcmp(logtable.squaretype, 'D') & (logtable.block == 3));
        
        matchmat.coord = nan (length(dRidx_b1), length(dRidx_b2)); 
        for iB = 1: 2    
            foo = eval(['logtable.path(dRidx_b' num2str(iB) ')']);
            for i = 1:length(foo)
                pathstr = foo{i};
                sep = strfind(foo{i}, ';'); 
                coordstr{iB, i} = pathstr(sep(end-1)+1 : sep(end)-1); 
            end 
            clear foo
        end  
        if ~isempty(dRidx_b3) 
            foo = logtable.path(dRidx_b3);
            for i = 1:length(foo)
                pathstr = foo{i};
                sep = strfind(foo{i}, ';'); 
                coordstr{3, i} = pathstr(sep(end-1)+1 : sep(end)-1); 
            end 
            clear foo        
        end 
        
        for idB1 = 1: length(dRidx_b1)
            for idB2 = 1: length(dRidx_b2)
                matchmat.coord(idB1, idB2) = strcmp(coordstr{1, idB1}, coordstr{2, idB2});
           end 
        end 
        
        % Find decision points that match in all 3 criteria 
%         matchmat.sum = matchmat.blackremains + matchmat.coord + matchmat.numsquaresopen;
        [i, j] = find(matchmat.coord == 1); 
        
       dt_b1 = nan(length(i), 1); dt_b2 = nan(length(i), 1); dt_b3 = nan(length(i), 1); 
       dpmatch(iW).matchingdp_b1 = zeros(size(dRidx_b1)); 
       dpmatch(iW).matchingdp_b2 = zeros(size(dRidx_b2)); 
       dpmatch(iW).matchingdp_b3 = zeros(size(dRidx_b3));       
       
        for iMatch = 1: length(i)           
            % converting to logtable row index to find the timestamps
            D_b1 = dRidx_b1(i(iMatch)); 
            D_b2 = dRidx_b2(j(iMatch));
            b3_match = [];
           
            decision_b1 = logtable.timestamp(D_b1); 
            postdecision_b1 = logtable.timestamp(D_b1+1); 
            dt_b1(iMatch) = postdecision_b1 - decision_b1;

            decision_b2 = logtable.timestamp(D_b2); 
            postdecision_b2 = logtable.timestamp(D_b2+1); 
            dt_b2(iMatch) = postdecision_b2 - decision_b2;
            
            if ~isempty(dRidx_b3) 
                % find b3 events that are at decision points that matched btw b1-b2   
                b3_match = find(double(cell2mat(cellfun(@(x) strcmp(coordstr{1, i(iMatch)}, x), {coordstr{3,:}}, 'UniformOutput', false)))); 
                if ~isempty(b3_match) 
                    D_b3 = dRidx_b3(b3_match); 
                    decision_b3 = logtable.timestamp(D_b3); 
                    postdecision_b3 = logtable.timestamp(D_b3+1); 
                    dt_b3(iMatch) = postdecision_b3 - decision_b3;
                    dpmatch(iW).matchingdp_b3(b3_match) = 1;   
                else 
                    dt_b3(iMatch) = NaN; 
                end         
            end 
            dpmatch(iW).matchingdp_b1(iMatch) = 1;  
            dpmatch(iW).matchingdp_b2(j(iMatch)) = 1;   
        end
        
        % mean decision time in the matching decision points
        data(iW,10) = mean(dt_b1);
        data(iW,11) = mean(dt_b2);
        data(iW,12) = mean(dt_b3);
         
        if mazeTidx_b3 == 0 
             if data(iW,10) <  data(iW,11)
                 data(iW,13) = 0;
             else
                 data(iW,13) = 1;
             end
        else
            if  data(iW,10) < data(iW,11) || data(iW,10) < data(iW,12)
                data(iW,13) = 0;
            else
                data(iW,13) = 1;
            end
        end        
        
        %% tidx per block
        data(iW,14) = mazeTidx_b1;
        data(iW,15) = mazeTidx_b2;
        data(iW,16) = mazeTidx_b3;
         %% save matching (b1-b2) decision points info 
        dpmatch(iW).mazename = worlds{iW};
        dpmatch(iW).tidx_b1 = mazeTidx_b1;
        dpmatch(iW).tidx_b2 = mazeTidx_b2;
        dpmatch(iW).tidx_b3 = mazeTidx_b3;
        %% Compare total moves taken in each block to the optimal number of moves it takes to reach the goal
        data(iW,17) = length(mazeSquaretype_b1) - maze_difficulty.optimal_steps(iW);
        data(iW,18) = length(mazeSquaretype_b2) - maze_difficulty.optimal_steps(iW);
        
        if mazeTidx_b3 == 0 
            data(iW,19) = 0;
        else
            data(iW,19) = length(mazeSquaretype_b3) - maze_difficulty.optimal_steps(iW);
        end
        
        % was there any improvement between blocks?           
        if mazeTidx_b3 == 0 
             if data(iW,17)  < data(iW,18)
                 data(iW,20) = 0;
             else
                 data(iW,20) = 1;
             end
        else
            if  data(iW,17) < data(iW,18) || data(iW,17) < data(iW,19)
                data(iW,20) = 0;
            else
                data(iW,20) = 1;
            end
        end        
               
        %Add optimal moves to data structure
        data(iW,21) = maze_difficulty.optimal_steps(iW);
        %% Compare maximum exploration path and optimal path to reveal goal square in each block to actual path taken
        % block 1
        squaresopen_b1 = logtable.numsquaresopen(find(strcmp(logtable.world, worlds(iW)) & (logtable.block == 1)));
        data(iW,22) = sum(squaresopen_b1 ~= 0) - maze_difficulty.optimal_path(iW);
        % block 2
        squaresopen_b2 = logtable.numsquaresopen(find(strcmp(logtable.world, worlds(iW)) & (logtable.block == 2)));
        data(iW,23) = sum(squaresopen_b2 ~= 0) - maze_difficulty.optimal_path(iW);
        % block 3
        squaresopen_b3 = logtable.numsquaresopen(find(strcmp(logtable.world, worlds(iW)) & (logtable.block == 3)));
         
        if mazeTidx_b3 == 0 
             data(iW,24) = 0;
        else
            data(iW,24) = sum(squaresopen_b3 ~= 0) - maze_difficulty.optimal_path(iW);
        end
        
        % b1-b2
        if data(iW,22)  < data(iW,23)
             data(iW,27) = 0;
        else
             data(iW,27) = 1;
        end     
         
        % b1-b3
        if mazeTidx_b3 == 0 
            data(iW,28) = 0;
        else 
            if  data(iW,22) < data(iW,24)
                data(iW,28) = 0;
            else
                data(iW,28) = 1;
            end
        end
         
        %Add optimal moves to data structure
        data(iW,25) = maze_difficulty.total_path(iW);
        data(iW,26) = maze_difficulty.optimal_path(iW);
        clear coordstr dRidx_*; 
    end 
    %% Finalize data table   
    data_reordered = cat(2, data(:, 1), data(:, 14:16), data(:, 2: 13),data(:, 17: 28)); 
    % convert behavioral data to table
    maze_nav_impv = array2table(data_reordered, 'VariableNames', {'world', 'TrialIdx_b1', 'TrialIdx_b2','TrialIdx_b3', ...
                                'TotalTime_b1', 'TotalTime_b2','TotalTime_b3', 'TotalTime_impv', ...
                                'nMoves_b1', 'nMoves_b2','nMoves_b3', 'nMoves_impv', ...
                                'DecisionTime_b1', 'DecisionTime_b2', 'DecisionTime_b3', 'DecisionTime_impv', ...
                                'OptimalStepOffset_b1','OptimalStepOffset_b2','OptimalStepOffset_b3','OptimalStepOffset_impv', 'nOptimalSteps',...
                                'OptimaPathOffset_b1','OptimaPathOffset_b2','OptimaPathOffset_b3', 'nTotalPath', 'nOptimalPath', 'OptimalPathOffset_impv12', 'OptimalPathOffset_impv13'});
    maze_nav_impv.world = worlds;
    % make sure maze names match in data and difficulty mat 
    if ~isequal(maze_difficulty.Var1, worlds) %words - is taken from logtable 
        error('Achtung: worlds order in data and difficulty file do not match!!'); 
    else 
        maze_nav_impv = [maze_nav_impv maze_difficulty.difficulty];
        maze_nav_impv.Properties.VariableNames([29]) = {'difficulty'};
        
    end 

    %% SAVE 
    save(fullfile(outdir, [subj, '_maze_nav_impv.mat']), 'maze_nav_impv'); 
    save(fullfile(outdir, [subj, '_matching_decpnts.mat']), 'dpmatch'); 
    
end 
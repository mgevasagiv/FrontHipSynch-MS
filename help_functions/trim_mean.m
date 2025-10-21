function [trimmed_mean, array, rm_bool] = trim_mean(array, by_std)
    m = nanmean(array);
    sd = std(array);
    rm_bool = (array > m + sd*by_std)| (array < m - sd*by_std); 
    array (rm_bool) = [];
    trimmed_mean = mean(array); 
end 

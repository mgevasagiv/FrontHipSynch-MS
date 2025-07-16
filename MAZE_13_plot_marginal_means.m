% Ensure blk and ripple_occur are categorical
data.blk = categorical(data.blk);
data.ripple_occur = categorical(data.ripple_occur);

% Create grid of all blk × ripple_occur combinations
blk_levels = categories(data.blk);
ripple_levels = categories(data.ripple_occur);

% Build new table for predictions
[blk_grid, ripple_grid] = ndgrid(blk_levels, ripple_levels);

newTbl = table;
newTbl.blk = categorical(blk_grid(:), blk_levels);
newTbl.ripple_occur = categorical(ripple_grid(:), ripple_levels);

% Grab a valid epair value from the original dataset
valid_epair = data.epair(1);  % or mode(data.epair) if you're being careful
valid_event = data.event(1);

% Make sure the newTbl.epair is the same type
newTbl.epair = repmat(valid_epair, height(newTbl), 1);

% Assign reference values for other fixed effects
% Choose the most common or representative value
newTbl.event = repmat(valid_event, height(newTbl), 1);

[mg_hat, mg_ci] = predict(lme_simple, newTbl);

blk_num = double(newTbl.blk);
ripple_labels = categories(newTbl.ripple_occur);

figure; hold on;
colors = lines(numel(ripple_labels));

for i = 1:numel(ripple_labels)
    idx = newTbl.ripple_occur == ripple_labels{i};
    
    % Jitter blk for visibility if necessary
    blk_vals = blk_num(idx);
    
    % Plot means with error bars
    errorbar(blk_vals, mg_hat(idx), ...
             mg_hat(idx) - mg_ci(idx,1), mg_ci(idx,2) - mg_hat(idx), ...
             '-o', 'Color', colors(i,:), 'DisplayName', ripple_labels{i});
end

xlabel('Block');
ylabel('Estimated mg');
legend('Location','best');
title('Estimated Marginal Means by Block and Ripple Occurrence');
axis([0 3 0 3])
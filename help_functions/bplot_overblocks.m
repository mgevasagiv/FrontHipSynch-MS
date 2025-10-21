function bplot_overblocks(data_b1, data_b2, txt_title, txt_ylabel)
    plot([data_b1 data_b2]', 'o-', 'LineWidth', 1.5, 'MarkerSize', 7);
    title(txt_title);
    xlim([0.5 2.5]);  xticks([1 2]); 
    xticklabels({'round1', 'round2'}); 
    ylabel(txt_ylabel);  set(gca,'fontsize', 18); 
end 
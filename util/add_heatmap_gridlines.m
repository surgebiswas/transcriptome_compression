function add_heatmap_gridlines( y )

set(gca, 'XTick', 0.5:size(y,2)+0.5);
set(gca, 'YTick', 0.5:size(y,1)+0.5);
set(gca, 'XTickLabel', []);
set(gca, 'YTickLabel', []);
set(gca, 'TickLength', [0 0]);
grid on;
set(gca, 'GridAlpha', 1);


end


function save_fig1_subheatmap( y, type, plottitle )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
tpmcmap = ((cbrewer('seq', 'YlOrRd', 64)));

figure
if strcmpi(type, 'gene')
    imagesc(sqrt(exp(y)), [0 3]); colormap(tpmcmap)
elseif strcmpi(type, 'geneset')
    imagesc(y, [-2 2]); colormap(parula)
end

axis equal;
axis tight;
daspect([1,2,1]);
add_heatmap_gridlines(y);
plotSave(plottitle);
close


end


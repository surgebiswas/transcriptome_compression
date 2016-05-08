clear;
rng('default')
cd ~/GitHub/data/transcriptome_compression/fig1/

path(genpath('~/GitHub/transcriptome_compression'), path)
bmus = [-2 0; 1 -1; 0 1; -3 2; 1 0];

b = [];
for i = 1 : size(bmus,1)
    b = [b; randn(5,2)*0.5 + repmat(bmus(i,:), 5, 1)];
end
b = b';

markers = randn(50,2);

y = standardize(markers*b);
ygs = standardize(markers*bmus');

ri=hclust(y);
y = y(ri,:);
ygs = ygs(ri,:);
markers = markers(ri,:); %markers(:,1) = -markers(:,1);


% Full gene expression matrix
save_fig1_subheatmap(y, 'gene', 'full_gene_heatmap.png');

% Gene set matrix
save_fig1_subheatmap(ygs, 'geneset', 'geneset_heatmap.png');


% Separated 
mreps = [5, 14]; 
save_fig1_subheatmap(y(:,1:mreps(1)-1), 'gene', 'sub_gene_heatmap_1.png');
save_fig1_subheatmap(y(:,mreps(1)+1:mreps(2)-1), 'gene', 'sub_gene_heatmap_2.png');
save_fig1_subheatmap(y(:,mreps(2)+1:end), 'gene', 'sub_gene_heatmap_3.png');

save_fig1_subheatmap(y(:,mreps(1)), 'gene', 'marker1_heatmap.png');
save_fig1_subheatmap(y(:,mreps(2)), 'gene', 'marker2_heatmap.png');

% Measurements -- decoding
subidx = 1:2:size(y,1);
save_fig1_subheatmap(-markers(subidx,1), 'gene', 'measured_marker1.png'); % negative because we're using col 5 in the full exp mat.
save_fig1_subheatmap(markers(subidx,2), 'gene', 'measured_marker2.png');


ygsh = standardize(markers*bmus' + 0.3*randn(size(ygs)));
save_fig1_subheatmap(ygsh(subidx,:), 'geneset', 'geneset_predicted_heatmap.png');

yh = standardize(markers*b + 0.3*randn(size(y)));
save_fig1_subheatmap(yh(subidx,:), 'gene', 'full_gene_predicted_heatmap.png');

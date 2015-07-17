function heatmap_raw_vs_reconstructed( Y, somp, model, organism, runclustering )
% Y = [samples x genes] standardized expression matrix, not in log-scale. 
% somp = result object from marker_OMP.
% model = tratrain model.
% organism = string specifying organism name.

DISPLAYRANGE = 2;
sY = standardize(Y);
lY = log(Y + 0.1);

if runclustering
    fprintf('Running clustering ... ');
    [ri, ci] = hclust(sY);
    save(sprintf('NCBI_SRA_%s_cluster_idx_for_raw_vs_reconstructed_heatmap.mat', organism), 'ri', 'ci');
    fprintf('Done.\n');
    return
else
    load(sprintf('NCBI_SRA_%s_cluster_idx_for_raw_vs_reconstructed_heatmap.mat', organism));
end


% Make percent variance explained vector. Sort by permuatation index in the clustered
% heatmap.
[~,~, ic] = intersect(somp.S(1:NUMMARKERS), ci, 'stable' );
p = diff([0, 1 - somp.punexp]);
pexp = zeros(1,length(ci));
for i = 1 : length(ic)
   pexp(ic(i)) = p(i);
end


% Reconstructed heatmap.
% Model reconstruction is in log scale
model.reconstruction = tradict(lY(:,somp.S), model);
sY_reconstructed = standardize(10.^(model.reconstruction) - 0.1); % undoes log10(x + 0.1), and then standardizes.

figure;
imagesc(sY_reconstructed(ri,ci), [-DISPLAYRANGE DISPLAYRANGE]);colormap(prgn);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'TickLength', [0 0]);
box on;
plotSave(sprintf('figures/heatmap_original_vs_reconstruction/%s_reconstructed_heatmap.png', organism));
close



% Original heatmap
figure;
imagesc(sY(ri,ci), [-DISPLAYRANGE DISPLAYRANGE]);colormap(prgn);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'TickLength', [0 0]);
box on;
plotSave(sprintf('figures/heatmap_original_vs_reconstruction/%s_original_heatmap.png', organism));
close

% Marker percent variance explained.
figure;
hold on;
qq = find(pexp);
for i = 1 : length(qq)
    plot([qq(i), qq(i)], [0 pexp(qq(i))], '-k', 'LineWidth', 3);
end
axis tight;
v = axis;
axis([0 length(pexp) 0 v(4)]);
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'TickLength', [0 0]);
box off;
plotSave(sprintf('figures/heatmap_original_vs_reconstruction/%s_marker_pexp_bars.png', organism));
close

% Colorbar
figure;
colorbar;
axis off
colormap(prgn);
caxis([-DISPLAYRANGE DISPLAYRANGE]);
plotSave('figures/heatmap_original_vs_reconstruction/colorbar.png');
close






end


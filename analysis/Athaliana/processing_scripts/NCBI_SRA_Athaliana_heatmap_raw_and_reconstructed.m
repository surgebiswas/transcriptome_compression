function NCBI_SRA_Athaliana_heatmap_raw_and_reconstructed( Y, somp, model )

COOVCUTOFF = 0;
NUMMARKERS = 100;

coov = std(Y')./mean(Y');
k = coov > COOVCUTOFF;
iset = find(k);

% Make sure that all of the markers are contained in the set of variably
% expressed genes.
assert( isempty(setdiff(somp.S(1:NUMMARKERS),iset)) );
sY = standardize(Y(iset,:)');

if false
    fprintf('Running clustering ... ');
    cg = clustergram(sY, 'Colormap', prgn, 'OptimalLeafOrder', false, 'DisplayRange', 3);
    save('NCBI_SRA_Athaliana_clustergram_for_raw_vs_reconstructed_heatmap.mat', 'cg');
    fprintf('Done.\n');
else
    load('NCBI_SRA_Athaliana_clustergram_for_raw_vs_reconstructed_heatmap.mat');
end

ri = str2double(get(cg, 'RowLabels'));
ci = str2double(get(cg, 'ColumnLabels'));

[~,~, ic] = intersect(somp.S(1:NUMMARKERS), ci, 'stable' );
p = diff([0, 1 - somp.punexp]);
pexp = zeros(1,length(ci));
for i = 1 : length(ic)
   pexp(ic(i)) = p(i);
end

% Reconstructed heatmap.
% Model reconstruction is in log scale
sY_reconstructed = standardize(10.^(model.reconstruction) - 0.1); % undoes log10(x + 0.1), and then standardizes.

figure;
imagesc(sY_reconstructed(ri,ci), [-2 2]);colormap(prgn);
%axis image
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'TickLength', [0 0]);
box on;
plotSave('figures/heatmap_original_vs_reconstruction/Athaliana_reconstructed_heatmap.png');
close



% Original heatmap
figure;
imagesc(sY(ri,ci), [-2 2]);colormap(prgn);
%axis image
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'TickLength', [0 0]);
box on;
plotSave('figures/heatmap_original_vs_reconstruction/Athaliana_original_heatmap.png');
close

% Marker percent variance explained.
figure;
hold on;
qq = find(pexp);
for i = 1 : length(qq)
    plot([qq(i), qq(i)], [0 pexp(qq(i))], '-k', 'LineWidth', 3);
   % plot(qq(i), pexp(qq(i)), 'ok', 'LineWidth', 3, 'MarkerFaceColor', 'k');
end
% set(h, 'FaceColor', 'k');
% set(h, 'EdgeColor', 'k');
axis tight;
v = axis;
axis([0 length(pexp) 0 v(4)]);
%axis off
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'TickLength', [0 0]);
box off;
plotSave('figures/heatmap_original_vs_reconstruction/marker_pexp_bars.png');
close



end


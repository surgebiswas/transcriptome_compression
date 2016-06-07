function  plot_perf_vs_num_markers

load('perf_vs_num_markers.mat')
load('perf_vs_num_markers_topN.mat');
load('perf_vs_num_markers_baseline_models.mat');
cm = lines;
sf = get_standard_figure_font_sizes;

% rng(1);
% sr_gene = log10(1:100)/3.1 + (0.1./(1:100)).*randn(1,100);  sr_gene(1) = 0.02; sr_gene = [0, sr_gene];
% lwa_gene = log10(1:100)/10 + (0.1./(1:100)).*randn(1,100) + 0.4; lwa_gene = [0, lwa_gene];
% 
% sr_genesets = log10(1:100)/4.3 + (0.1./(1:100)).*randn(1,100) + 0.45;  sr_genesets(1) = 0.05; sr_genesets = [0, sr_genesets];
% lwa_genesets = log10(1:100)/10 + (0.1./(1:100)).*randn(1,100) + 0.65; lwa_genesets = [0, lwa_genesets];



clf
subplot(1,2,1)
hold on
l(1) = plot(0:100, [0 pcc_genes], '-k', 'LineWidth', 3);
l(2) = plot(0:100, [0 pcc_genes_topN], '-', 'LineWidth', 3, 'Color', [0.5 0.5 0.5]);
l(3) = plot(0:100, sr_gene, '-', 'LineWidth', 3, 'Color', cm(1,:));
l(4) = plot(0:100, lwa_gene, '-', 'LineWidth', 3, 'Color', cm(2,:));
axis square
axis([0 100 0 1])
set(gca, 'YTick', 0:0.1:1);
set(gca, 'XTick', [0 25 50 75 100]);
set(gca, 'FontSize', sf.axis_tick_labels-4);
xlabel('Num. markers', 'FontSize', sf.axis_labels-4);
ylabel('Avg. PCC', 'FontSize', sf.axis_labels-4);
sub_pos = get(gca,'position'); % get subplot axis position
set(gca,'position',sub_pos.*[1 1 1.18 1.18]) % stretch its width and height
title('genes', 'FontSize', sf.axis_labels-4);




subplot(1,2,2)
hold on
l(1) = plot(0:100, [0 pcc_genesets], '-k', 'LineWidth', 3);
l(2) = plot(0:100, [0 pcc_genesets_topN], '-', 'LineWidth', 3, 'Color', [0.5 0.5 0.5]);
l(3) = plot(0:100, sr_genesets, '-', 'LineWidth', 3, 'Color', cm(1,:));
l(4) = plot(0:100, lwa_genesets, '-', 'LineWidth', 3, 'Color', cm(2,:));
legend(l, 'Tradict', 'Tradict Shallow-Seq', 'SR', 'LWA', 'Location', 'SouthEast');
axis square
axis([0 100 0 1])
set(gca, 'XTick', [0 25 50 75 100]);
set(gca, 'FontSize', sf.axis_tick_labels-4);
xlabel('Num. markers', 'FontSize', sf.axis_labels-4);

set(gca, 'YTick', []);
sub_pos = get(gca,'position'); % get subplot axis position
set(gca,'position',sub_pos.*[1 1 1.18 1.18]) % stretch its width and height
title('tr. programs', 'FontSize', sf.axis_labels-4);

plotSave('figures/perf_v_num_markers.png');
close

end


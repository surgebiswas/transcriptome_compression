clear;
rng('default');
cd('/Users/sbiswas/GitHub/data/transcriptome_compression/fig2');

load('/Users/sbiswas/GitHub/data/transcriptome_compression/Mmusculus/NCBI_SRA_Mmusculus_PCA_pexp_vs_eigengene_params.mat');
mpexp = pexp;
mpexpn = pexp_null;

load('/Users/sbiswas/GitHub/data/transcriptome_compression/Athaliana/NCBI_SRA_Athaliana_PCA_pexp_vs_eigengene_params.mat');
apexp = pexp;
apexpn = pexp_null;

pexp = [apexp(1:500), mpexp(1:500)];
pexpn = [apexpn(1:500), mpexpn(1:500)];


sf = get_standard_figure_font_sizes;
figure;
hold on
s = semilogx(cumsum(pexp), '-', 'LineWidth', 3);
semilogx(cumsum(pexpn(:,1)), '--', 'LineWidth', 3, 'Color', s(1).Color);
semilogx(cumsum(pexpn(:,2)), '--', 'LineWidth', 3, 'Color', s(2).Color)
set(gca, 'XScale', 'log');
set(gca, 'XTick', [0 1 10 100 500]);
set(gca, 'XTickLabel', [0 1 10 100 500]);
axis([0 500 0 100]);
grid on
axis square
box on
legend(s, 'A. thaliana', 'M. musculus', 'Location', 'NorthWest');
xlabel('Num. components', 'FontSize', sf.axis_labels);
ylabel('Percent variance explained', 'FontSize', sf.axis_labels);
set(gca, 'FontSize', sf.axis_tick_labels);
set(gca, 'GridAlpha', 0.7);
plotSave('pexp_v_num_components.png');
close;



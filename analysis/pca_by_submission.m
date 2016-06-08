function pca_by_submission( s, qt, pexp )
rng('default')
sf = get_standard_figure_font_sizes;

%s = lY*coef(:,1:2);

figure;
hold on
usub = unique(qt.Submission);
covs = cell(length(usub),1);
covnorms = zeros(length(usub),1);
for i = 1 : length(usub)
    idx = strcmpi(qt.Submission, usub{i});
    plot(s(idx,1), s(idx,2), '.', 'Color', rand(1,3), 'MarkerSize', 6);
    
    covs{i} = cov(s(idx,:));
    covnorms(i) = norm(covs{i}, 'fro');
end

box on;
axis square;
axis tight;
buffer_axis;
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'TickLength', [0 0 ]);

xlabel(sprintf('PC1 (%0.1f%%)', pexp(1)), 'FontSize', sf.axis_labels);
ylabel(sprintf('PC2 (%0.1f%%)', pexp(2)), 'FontSize', sf.axis_labels);




mc = median(covnorms);
[~, mcovidx] = min( abs(covnorms - mc) );
mcov = diag(diag(covs{mcovidx}));


% Calculate the location of the ellipse. Put it in the bottom left.
v = axis;
normloc = 0.1;
xl = v(2) - v(1);
yl = v(4) - v(3);

xloc = v(1) + (1 - normloc)*xl;
yloc = v(3) + normloc*yl;


h = plot_gaussian_ellipsoid([xloc yloc], mcov);
set(h, 'Color', 'k');
set(h, 'LineWidth', 3);
plotSave('figures/pca/pca_by_submission.png');
close

 



end


function NCBI_SRA_Mmusculus_pexp_vs_components( sY )
% sY = lY (don't standardize lY)
    
    sf = get_standard_figure_font_sizes;
    NUMCOMPONENTS = 500;
    LEGENDFONTSIZE = 10;

    sYnull = sY;
    for i = 1 : size(sY,2)
        sYnull(:,i) = sY(randperm(size(sY,1)),i);
    end
    
    if false
        disp('Computing for full model');
        [coef, ~, ~, ~, pexp] = pca(sY, 'NumComponents', NUMCOMPONENTS);
        
        disp('Computing for null model');
        [~, ~, ~,~, pexp_null] = pca(sYnull, 'NumComponents', NUMCOMPONENTS);
        save('NCBI_SRA_Mmusculus_PCA_pexp_vs_eigengene_params.mat', 'coef', 'pexp', 'pexp_null');
    else
        load('NCBI_SRA_Mmusculus_PCA_pexp_vs_eigengene_params.mat');
    end 
    cp = cumsum(pexp);
    cp_null = cumsum(pexp_null);

    
    % Copied from A. thaliana code.
    figure;
    h(1) = semilogx(cp(1:NUMCOMPONENTS), '-k', 'LineWidth', 3);
    hold on
    h(2) = semilogx(cp_null(1:NUMCOMPONENTS), '--k', 'LineWidth', 3);
    plot([100 100], [0 cp(100)], '-r', 'LineWidth', 2);
    hold on
    plot([1 100], [cp(100), cp(100)], '-r', 'LineWidth', 2);
    l = legend(h, 'PCA - Original data', 'PCA - Permuted data', 'Location', 'NorthWest');
    set(l, 'FontSize', sf.axis_tick_labels);
    set(gca, 'YTick', 0:10:100);
    axis square
    set(gca, 'FontSize', sf.axis_tick_labels);
    axis tight;
    v = axis;
    axis([v(1) v(2) 0 100]);
    xlabel('Principal component', 'FontSize', sf.axis_labels);
    ylabel('Percent variance explained', 'FontSize', sf.axis_labels);

    
    plotSave('figures/pca/NCBI_SRA_Mmusculus_Pexp_vs_eigengene_full.png');
    iminvert('figures/pca/NCBI_SRA_Mmusculus_Pexp_vs_eigengene_full.png');
    
    
    
    if true
        load('NCBI_SRA_Mmusculus_marker_OMP_decomposition_punexp_0.00_maxfeats_500.mat');
        hold on
        h(3) = plot(100 - 100*somp.punexp, '-g', 'LineWidth', 3);
        plot([100 100], [0 cp(100)], '-r', 'LineWidth', 2); % Re do vertical red line.
        plot([1 100], 100 - 100*somp.punexp(100)*[1 1], '-r', 'LineWidth', 2);
        l = legend(h, 'PCA - Original data', 'PCA - permuted data', 'Tradict - Original data', 'Location', 'NorthWest');
        set(l, 'FontSize', sf.axis_tick_labels);
        xlabel('Number of features', 'FontSize', sf.axis_labels);
        
        plotSave('figures/pca/NCBI_SRA_Mmusculus_Pexp_vs_eigengene_full_with_MOMP.png');
        iminvert('figures/pca/NCBI_SRA_Mmusculus_Pexp_vs_eigengene_full_with_MOMP.png');
        close
        
       
    end

end


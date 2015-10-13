function [ coef, pexp ] = pexp_vs_components( sY, organism )

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
        s = sY*coef;
        save(sprintf('NCBI_SRA_%s_PCA_pexp_vs_eigengene_params.mat', organism), 'coef', 'pexp', 'pexp_null');
        return
    else
        load(sprintf('NCBI_SRA_%s_PCA_pexp_vs_eigengene_params.mat', organism));
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
    l = legend(h, 'PCA - Original', 'PCA - Permuted', 'Location', 'NorthWest');
    set(l, 'FontSize', sf.axis_tick_labels);
    set(gca, 'YTick', 0:10:100);
    axis square
    set(gca, 'FontSize', sf.axis_tick_labels);
    axis tight;
    v = axis;
    axis([v(1) v(2) 0 100]);
    xlabel('Principal component', 'FontSize', sf.axis_labels);
    ylabel('Percent variance explained', 'FontSize', sf.axis_labels);

    
    plotSave(sprintf('figures/pca/NCBI_SRA_%s_Pexp_vs_eigengene_full.png', organism));
    iminvert(sprintf('figures/pca/NCBI_SRA_%s_Pexp_vs_eigengene_full.png', organism));
    
    
    
    if exist(sprintf('NCBI_SRA_%s_marker_OMP_decomposition_punexp_0.00_maxfeats_500.mat', organism), 'file')
        
        load(sprintf('NCBI_SRA_%s_marker_OMP_decomposition_punexp_0.00_maxfeats_500.mat', organism));
        hold on
        h(3) = plot(100 - 100*somp.punexp, '-g', 'LineWidth', 3);
        plot([100 100], [0 cp(100)], '-r', 'LineWidth', 2); % Re do vertical red line.
        plot([1 100], 100 - 100*somp.punexp(100)*[1 1], '-r', 'LineWidth', 2);
        l = legend(h, 'PCA - Original', 'PCA - Permuted', 'Tradict - Original', 'Location', 'NorthWest');
        set(l, 'FontSize', sf.axis_tick_labels);
        xlabel('Number of features', 'FontSize', sf.axis_labels);
        
        plotSave(sprintf('figures/pca/NCBI_SRA_%s_Pexp_vs_eigengene_full_with_MOMP.png', organism));
        iminvert(sprintf('figures/pca/NCBI_SRA_%s_Pexp_vs_eigengene_full_with_MOMP.png', organism));
        close
        
       
    end

    

end


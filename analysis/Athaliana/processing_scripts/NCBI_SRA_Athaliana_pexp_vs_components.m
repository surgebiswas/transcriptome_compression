function [ coef, pexp ] = NCBI_SRA_Athaliana_pexp_vs_components( sY )

%sY = An expression table. samples x genes.
    sf = get_standard_figure_font_sizes;
    
    NUMCOMPONENTS = 500;

    sYnull = sY;
    for i = 1 : size(sY,2)
        sYnull(:,i) = sY(randperm(size(sY,1)),i);
    end
    
    if false
        [coef, ~, ~, ~, pexp] = pca(sY, 'NumComponents', NUMCOMPONENTS);
        [~, ~, ~,~, pexp_null] = pca(sYnull, 'NumComponents', NUMCOMPONENTS);
        save('PCA_pexp_vs_eigengene_params.mat', 'coef', 'pexp', 'pexp_null');
    else
        load('PCA_pexp_vs_eigengene_params.mat');
    end 
    cp = cumsum(pexp);
    cp_null = cumsum(pexp_null);

    
    
    figure;
    h(1) = semilogx(cp(1:NUMCOMPONENTS), '-k', 'LineWidth', 3);
    hold on
    h(2) = semilogx(cp_null(1:NUMCOMPONENTS), '--k', 'LineWidth', 3);
    plot([100 100], [0 cp(100)], '-r', 'LineWidth', 2);
    hold on
    plot([1 100], [cp(100), cp(100)], '-r', 'LineWidth', 2);
    legend(h, 'Original data', 'Null model', 'Location', 'NorthWest');
    axis square
    set(gca, 'FontSize', sf.axis_tick_labels);
    axis tight;
    xlabel('Principal component (eigengene)', 'FontSize', sf.axis_labels);
    ylabel('Percent variance explained', 'FontSize', sf.axis_labels);
    grid on
    
    plotSave('figures/Pexp_vs_eigengene_full.png');
    iminvert('figures/Pexp_vs_eigengene_full.png');
    
    
    % Add the marker OMP decomposition results
    if true
        load('NCBI_SRA_Athaliana_marker_OMP_decomposition_punexp_0.00_maxfeats_500.mat');
        hold on
        h(3) = plot(100 - 100*somp.punexp, '-c', 'LineWidth', 3);
        plot([100 100], [0 cp(100)], '-r', 'LineWidth', 2); % Re do vertical red line.
        l = legend(h, 'PCA', 'Null model', 'Tradict', 'Location', 'NorthWest');
        set(l, 'FontSize', 12);
        xlabel('Number of features', 'FontSize', sf.axis_labels);
        grid on
        
        plotSave('figures/Pexp_vs_eigengene_full_with_MOMP.png');
        iminvert('figures/Pexp_vs_eigengene_full_with_MOMP.png');
        close
        
        
        % Marker OMP results only. This figure is being used for the talk
        % given to David Sainsbury.
%         figure
%         h2(1) = semilogx(100 - 100*somp_full.punexp, '-c', 'LineWidth', 3);
%         hold on;
%         h2(2) = semilogx(cp_null(1:length(somp_full.punexp)), '--k', 'LineWidth', 3);
%         legend(h2, 'Tradict', 'Null model', 'Location', 'NorthWest');
%         grid on;
%         axis square
%         set(gca, 'FontSize', 12);
%         axis tight;
%         xlabel('Number of marker genes selected', 'FontSize', 14);
%         ylabel('Percent variance explained', 'FontSize', 14);
%         plotSave('figures/Tradict_pexp_vs_num_genes_no_PCA.png');
%         iminvert('figures/Tradict_pexp_vs_num_genes_no_PCA.png');
%         close
%         
        
    end

    
    close

end


function [ coef, pexp ] = NCBI_SRA_Athaliana_pexp_vs_components( sY )

%sY = An expression table. samples x genes.
    
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
    set(gca, 'FontSize', 12);
    axis tight;
    xlabel('Principal component (eigengene)', 'FontSize', 14);
    ylabel('Percent variance explained', 'FontSize', 14);
    plotSave('figures/Pexp_vs_eigengene_full.png');
    iminvert('figures/Pexp_vs_eigengene_full.png');
    
    
    % Add the marker OMP decomposition results
    if true
        load('NCBI_SRA_Athaliana_context_specific_performance_MOMP.mat')
        hold on
        h(3) = plot(100 - 100*somp_full.punexp, '-c', 'LineWidth', 3);
        legend(h, 'Original data', 'Null model', 'Marker OMP', 'Location', 'NorthWest');
        xlabel('Number of features', 'FontSize', 14);
        grid on
        plotSave('figures/Pexp_vs_eigengene_full_with_MOMP.png');
        iminvert('figures/Pexp_vs_eigengene_full_with_MOMP.png');
        close
        
        
        % Marker OMP results only. This figure is being used for the talk
        % given to David Sainsbury.
        figure
        h2(1) = semilogx(100 - 100*somp_full.punexp, '-c', 'LineWidth', 3);
        hold on;
        h2(2) = semilogx(cp_null(1:length(somp_full.punexp)), '--k', 'LineWidth', 3);
        legend(h2, 'Tradict', 'Null model', 'Location', 'NorthWest');
        grid on;
        axis square
        set(gca, 'FontSize', 12);
        axis tight;
        xlabel('Number of marker genes selected', 'FontSize', 14);
        ylabel('Percent variance explained', 'FontSize', 14);
        plotSave('figures/Tradict_pexp_vs_num_genes_no_PCA.png');
        iminvert('figures/Tradict_pexp_vs_num_genes_no_PCA.png');
        close
        
        
    end

    
    close

end


function NCBI_SRA_Mmusculus_pexp_vs_components( sY )
% sY = lY (don't standardize lY)

    NUMCOMPONENTS = 500;

    sYnull = sY;
    for i = 1 : size(sY,2)
        sYnull(:,i) = sY(randperm(size(sY,1)),i);
    end
    
    if true
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
    plotSave('figures/pca/Pexp_vs_eigengene_full.png');
    iminvert('figures/pca/Pexp_vs_eigengene_full.png');


end


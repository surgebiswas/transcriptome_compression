function NCBI_SRA_Mmusculus_plot_PCA( lY, coef, qt, pexp )
    
    sf = get_standard_figure_font_sizes;
    MSIZE = 15;
    UNANNOTSCALE = 5;
    LEGENDFONTSIZE = 10;
    SAVEFORBLACKBACKGROUND = true;


    LIGHTEN = 50;
    
    % PC scores
    s = lY*coef(:,1:3);
    
    
    
    
    queries = { {'hematopoetic', 'lymphatic'},  {'stem_cell'}, {'reproductive'}, {'breast_cancer'} ...
                {'embryonic', 'developing_liver', 'developing_cardiovascular', 'developing_respiratory', 'developing_pancreas', 'developing_urinary'}, ...
                {'connective', 'epithelium', 'integumentary', 'mammary'}, ...
                {'digestive', 'pancreas', 'cardiovascular', 'respiratory', 'urinary', 'thymus'}, ...
                {'skeletal', 'muscular'}, ...
                {'liver'}, ...
                {'nervous', 'retina'}, {'developing_nervous'}, ...
                {'adrenal', 'placenta', 'unannotated', 'NA'} };
                
    labels = {'hematopoetic/lymphatic', 'stem cell', 'reproductive', 'breast cancer', 'embryonic', ...
              'connective/epithelium/skin', ...
              'viscera', ...
              'musculoskeletal', 'liver', 'nervous', 'developing nervous', ...
              'annotation pending'};
    
    colors = [ 26, 255, 83; 220 20 60; 250 128 114; 255 0 255; 128 0 0; 
               255 215 0; 255 140 0; 189 183 107; 107 142 35; 69 40 255; 129 60 255;
               178 178 178]/255;
               
    %[ cbrewer('qual', 'Paired', length(queries)-1, 'cubic'); [0.7 0.7 0.7] ];
    
    figure;
    hold on;
    for i = 1 : length(queries)
        k = steq(qt.class, queries{i});
        
        if strcmpi(labels{i}, 'annotation pending')
            msize = MSIZE/UNANNOTSCALE;
        else
            msize = MSIZE;
        end
        
        if SAVEFORBLACKBACKGROUND
            h(i) = scatter3(s(k,1), s(k,2), s(k,3), msize, 'MarkerFaceColor', 1 - colors(i,:), ...
                'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
        else
            h(i) = scatter3(s(k,1), s(k,2), s(k,3), msize, 'MarkerFaceColor', colors(i,:), ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        end
    end
    l = legend(h, labels, 'Location','NorthEastOutside');
    box on;
    axis square;
    xlabel(sprintf('PC1 (%0.1f%%)', pexp(1)), 'FontSize', sf.axis_labels);
    ylabel(sprintf('PC2 (%0.1f%%)', pexp(2)), 'FontSize', sf.axis_labels);
    zlabel(sprintf('PC3 (%0.1f%%)', pexp(3)), 'FontSize', sf.axis_labels);
    set(l, 'Box', 'off');
    set(l, 'FontSize', LEGENDFONTSIZE);
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca, 'ZTick', []);
    set(gca, 'TickLength', [0 0]);
    
    
    viewpoints = [0 90; 90 0; 0 0];
    for i = 1 : size(viewpoints,1)
        view(viewpoints(i,:));
        axis tight;
        buffer_axis;
        if SAVEFORBLACKBACKGROUND
            plotSave(sprintf('figures/pca/NCBI_SRA_Mmusculus_PCA_az_%0.0f_el_%0.0f_for_black.png', viewpoints(i,1), viewpoints(i,2)));
            iminvert(sprintf('figures/pca/NCBI_SRA_Mmusculus_PCA_az_%0.0f_el_%0.0f_for_black.png', viewpoints(i,1), viewpoints(i,2)));
        else  
            plotSave(sprintf('figures/pca/NCBI_SRA_Mmusculus_PCA_az_%0.0f_el_%0.0f.png', viewpoints(i,1), viewpoints(i,2)));
        end
    end
    
    close


end


function NCBI_SRA_Mmusculus_plot_PCA( lY, coef, qt, pexp )
    
    sf = get_standard_figure_font_sizes;


    % PC scores
    s = lY*coef(:,1:3);
    
    
    MSIZE = 30;
    
    queries = { {'unannotated', 'NA'}, {'immune'}, {'lymphatic'}, {'hematopoetic'}, {'developing_central_nervous'}, ...
        {'central_nervous'}, {'embryo'}, {'liver'}, {'digestive'}, {'stem_cell'}, {'integumentary'}, ...
        {'reproductive'}, {'cardiovascular'}, {'respiratory'}, {'muscular'}, ...
        {'pancreas'}, {'retina'}, {'skeletal', 'urinary'} };
    labels = {'annotation pending', 'immune', 'lymphatic', 'hematopoetic', 'developing central nervous', ...
        'mature central nervous', 'embryo', 'liver', 'digestive', 'stem_cell', 'integumentary', ...
        'reproductive', 'cardiovascular', 'respiratory', 'muscular', 'pancreas', ...
        'retina', 'other'};
    
    colors = [ [0.7 0.7 0.7]; cbrewer('qual', 'Paired', length(queries)-1, 'cubic') ];
    
    
    figure;
    hold on;
    for i = 1 : length(queries)
        k = steq(qt.class, queries{i});
        
        if strcmpi(labels{i}, 'annotation pending')
            msize = MSIZE/2;
        else
            msize = MSIZE;
        end
        
        h(i) = scatter3(s(k,1), s(k,2), s(k,3), msize, 'MarkerFaceColor', colors(i,:), ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    end



end


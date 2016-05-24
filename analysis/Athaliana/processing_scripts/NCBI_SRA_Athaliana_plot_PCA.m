function NCBI_SRA_Athaliana_plot_PCA( lY, coef, qt, pexp )

sf = get_standard_figure_font_sizes;
MSIZE = 15;
UNANNOTSCALE = 5;
LEGENDFONTSIZE = 10;
SAVEFORBLACKBACKGROUND = false;

md = dataset('XLSfile', 'NCBI_SRA_Athaliana_run_metadata.xlsx', 'Sheet', 1, 'ReadObsNames', true, 'ReadVarNames', true);


% Link metadata with query table.
omd = get(md, 'ObsNames');
oqt = get(qt, 'ObsNames');
qt.tissue = cell(size(qt,1),1);
for j = 1 : length(qt.tissue);
    qind = find(strcmpi(oqt{j}, omd));
    if isempty(qind)
        qt.tissue{j} = 'unannotated';
    else
        assert(length(qind) == 1)
        qt.tissue{j} = md.tissue{qind};
    end
end




% PC scores
s = lY*coef(:,1:3);


    MINMSIZE = MSIZE;
    MAXMSIZE = MSIZE;

    os = 3; % artifact -- clean up later.
    maxos = max(s(:,os));
    minos = min(s(:,os));

    embryonic = {'developing seed', 'endosperm'};
    flower = {'flower', 'floral bud', 'carpel'};
    aground = {'leaves', 'shoot'};
    bground = {'root'};
    seedling = {'seedling'};
    unannotated = {'unannotated'};

    tissues = {embryonic, flower, aground, bground, seedling, unannotated};
    colors = {[210,22,30]/256, [218,165,32]/256, [15 243 61]/256, ...
        [243, 236, 15]/256, [15 150 159]/256, [0.7 0.7 0.7]};

    figure;
    hold on

    for i = 1 : length(tissues)
        k = steq(qt.tissue, tissues{i});
        norm_msizes =  (s(k,os) - minos)./(maxos - minos);
        msizes = (MAXMSIZE - MINMSIZE)*norm_msizes + MINMSIZE;
        
        if strcmpi(tissues{i}{1}, 'unannotated')
            msizes = msizes/UNANNOTSCALE;
        end
        
        if SAVEFORBLACKBACKGROUND
            h(i) = scatter3(s(k,1), s(k,2), s(k,3), msizes, 'MarkerFaceColor', 1 - colors{i}, ...
                'MarkerEdgeColor', 'w', 'LineWidth', 0.5);
        else
            h(i) = scatter3(s(k,1), s(k,2), s(k,3), msizes, 'MarkerFaceColor', colors{i}, ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        end


    end

    %grid on
    axis square
    box on
    xlabel(sprintf('PC1 (%0.1f%%)', pexp(1)), 'FontSize', sf.axis_labels);
    ylabel(sprintf('PC2 (%0.1f%%)', pexp(2)), 'FontSize', sf.axis_labels);
    zlabel(sprintf('PC3 (%0.1f%%)', pexp(3)), 'FontSize', sf.axis_labels);
    l = legend(h, {'seed/endosperm', 'flower/floral bud/carpel', 'leaves/shoot', 'root', 'seedling', 'annotation pending'}, 'Location', 'NorthEastOutside');
    set(l, 'Box', 'off');
    set(l, 'FontSize', LEGENDFONTSIZE);
    set(gca, 'FontSize', sf.axis_tick_labels);
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca, 'ZTick', []);
    set(gca, 'TickLength', [0 0]);
    
    
    % 1 vs 2, 2 vs 3, 1 vs 3, mix1
    viewpoints = [0 90; 90 0; 0 0];
    for i = 1 : size(viewpoints,1)
        view(viewpoints(i,:));
        axis tight;
        buffer_axis;
        
        if SAVEFORBLACKBACKGROUND
            plotSave(sprintf('figures/pca/NCBI_SRA_Athaliana_PCA_az_%0.0f_el_%0.0f_for_black.png', viewpoints(i,1), viewpoints(i,2)));
            iminvert(sprintf('figures/pca/NCBI_SRA_Athaliana_PCA_az_%0.0f_el_%0.0f_for_black.png', viewpoints(i,1), viewpoints(i,2)));
        else
            plotSave(sprintf('figures/pca/NCBI_SRA_Athaliana_PCA_az_%0.0f_el_%0.0f.png', viewpoints(i,1), viewpoints(i,2)));
        end
        
    end
    
    close
    
        
            


end


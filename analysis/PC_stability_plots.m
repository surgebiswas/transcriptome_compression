function PC_stability_plots( coef, pexp, lY, qt, varargin )
% Evolution of the gene expression space over time.
    sf = get_standard_figure_font_sizes;
    makemovie = setParam(varargin, 'makemovie', false);
    makeheatmaps = setParam(varargin, 'makeheatmaps', false);
    makersqplots = setParam(varargin, 'makersqplots', false);

    
    dn = datenum(qt.release_date);
    dnu = unique(dn); % sorted earliest to latest.
    md = min(dnu);
    xt = max(dnu) - min(dnu) : round(-365/1) : 0;
    
    s = lY*coef(:,1:10); % project onto the first few principal components.
    
    % Movie of PCA over time.
    if makemovie
        MAXCOMP = 3;
        cm = cbrewer('qual', 'Set1', length(dnu), 'cubic');
        figure;
        F(1) = struct('cdata',[],'colormap',[]); 
        set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); 
        
        for j = 1 : length(dnu)
            k = dn == dnu(j);
            
            sb = 1;
            for l = 1 : 2
                for m = l + 1 : 3
                    subplot(1,3,sb)
                    hold on
                    plot(s(k,l), s(k,m), '.', 'Color', cm(j,:));
                    drawnow
                    axis([min(s(:,l)), max(s(:,l)), min(s(:,m)), max(s(:,m))]);
                    buffer_axis;
                    box on
                    axis square;
                    xlabel(sprintf('PC%0.0f (%0.1f%%)', l, pexp(l)), 'FontSize', 14)
                    ylabel(sprintf('PC%0.0f (%0.1f%%)', m, pexp(m)), 'FontSize', 14)
                    set(gca, 'FontSize', 12);
                    
                    if sb == 2; title(datestr(dnu(j)), 'FontSize', 16); end
                    
                    
                    F(end+1) = getframe(gcf);
                    sb = sb+1;
                end
            end
        end
        F(1) = [];
        
        writerObj = VideoWriter('figures/PC_stability/PCA_over_time.avi');
        writerObj.FrameRate = 8;
        open(writerObj)

        writeVideo(writerObj, F);

        close(writerObj);
        close;
    end
    
    % Compute temporal stability maps
    mfile = 'figures/PC_stability/PC_stability_maps.mat';
    if ~exist(mfile, 'file');
        for j = 1 : 10; %length(pm)
            fprintf('Computing PC stability map %0.0f\n', j);
            [pm{j}, cmax(j)] = pc_temporal_density_map(s(:,j), dn, dnu);
        end
        save(mfile, 'pm', 'cmax');
    else
        load(mfile);
    end
    
    
    % Plot heatmaps
    if makeheatmaps
        cm = cbrewer('seq', 'YlOrRd', length(dnu), 'cubic');
        cmaxmax = 0.9*max(cmax);
        
        % colorbar
        cc = colorbar;
        axis off;
        colormap(cm);
        set(cc, 'TickLength', [0 0]);
        set(gca, 'CLim', [0 cmaxmax]);
        yt = linspace(0, cmaxmax, 4);
        set(cc, 'YTick', yt);
        set(cc, 'YTickLabel', round((10.^yt)*10)/10 - 1);
        set(gca, 'FontSize', 22)
        plotSave('figures/PC_stability/PC_temporal_density_colorbar.png')
        close
        
        
        for j = 1 : length(pm)
            figure;
            imagesc(flipud(pm{j}), [0, cmaxmax]); colormap(cm);
            
            axis image
            set(gca, 'TickLength', [0 0]);
            
            
            set(gca, 'XTick', fliplr(xt));
            set(gca, 'XTickLabel', flipud(cellstr(datestr(xt + min(dnu)))) ) ;
            rotateXLabels_imagesc(gca, 36);
            
            ytl = [min(s(:,j)), max(s(:,j))];
            v = axis;
            yt = [v(3), v(4)];
            set(gca, 'YTick', []);
            %set(gca, 'YTick', yt);
            %set(gca, 'YTickLabel', round(100*fliplr(ytl))/100);
            
            yl = ylabel(sprintf('PC%0.0f (%0.1f%%)', j, pexp(j)), 'FontSize', sf.axis_labels - 6);
            set(yl, 'Units', 'Normalized', 'Position', [-0.03, 0.5, 0]);
            
            set(gca, 'FontSize', sf.axis_tick_labels - 6);
            
            plotSave(sprintf('figures/PC_stability/PC%0.0f_temporal_density_plot.png', j));
            close
        end
        
       
        
        
        
         ims = [];   
         MAXCOMP = 3;
         for j = 1 : MAXCOMP
             im = imread(sprintf('figures/PC_stability/PC%0.0f_temporal_density_plot.png', j));
             if j == MAXCOMP
                im([1:370, 1170:end],:, :) = [];
             else
                 im([1:370, 880:end],:, :) = [];
             end
             ims = cat(1, ims, im);
         end
         imwrite(ims, sprintf('figures/PC_stability/PC%0.0f-%0.0f_temporal_density_plot_aggregate.png', 1, MAXCOMP));
    end
    
    % Plot R^2 to final over time.
    if makersqplots
        figure;
        hold on
        
        cm = cbrewer('qual', 'Set1', length(pm), 'cubic');
        for j = 1 : length(pm)
            rs = corr(pm{j}, pm{j}(:,end));
            plot(rs, '-', 'Color', cm(j,:), 'LineWidth', 2);
        end
        axis square;
        axis tight;
        box on;
        set(gca, 'XTick', fliplr(xt));
        set(gca, 'XTickLabel', flipud(cellstr(datestr(xt + min(dnu)))) ) ;
        rotateXLabels(gca, 36);
        set(gca, 'FontSize', sf.axis_tick_labels - 6);

        ylabel('Pearson correlation', 'FontSize', sf.axis_labels);
        plotSave('figures/PC_stability/PC_density_Rsq_vs_time.png');
        close
    end
      
    
    
    function [pm,cmax] = pc_temporal_density_map(si, dn, dnu)
        
        mind = min(dnu);
        maxd = max(dnu);
        daterange = maxd-mind;
        NQUERY = round(daterange/4);
        
        
        buffer = (max(si) - min(si))*0.05;
        query = linspace(min(si) - buffer, max(si) + buffer, NQUERY)';
        
        pm = zeros(NQUERY, daterange);
        for i = 1 : length(dnu)
            q = dn <= dnu(i);
            startIdx = dnu(i) - mind + 1;
            
            
            
            ff = ksdensity(si(q), query);
            
            
            fft = ff*sum(q)/sum(ff);
            pm(:,startIdx:end) = repmat(log10(fft+1), 1, size(pm(:,startIdx:end),2));
            cmax = max(max(pm));
        end
        
        
        
        
        
    end


end


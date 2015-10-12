function plot_quality_filter_report( reportFile, qfparams, organism )

rep = load(reportFile);

mr_cutoff = qfparams.MRTHRESH; %0.7;
md_cutoff = qfparams.RCTHRESH; %4e6;
corr_cutoff = qfparams.CORRCUTOFF; % 0.7
nz_cutoff = qfparams.NZCUTOFF; % 0.35
     


% Plot the mapping rate vs mapped depth plot
plot_mr_vs_md(rep.d,rep.mr);

% Plot average correlation to other samples vs. Prop non-zero expressed
% genes.
plot_ac_vs_pnz(rep.q, rep.z);

% Plot coefficient of variation density;
plot_coov_density(rep.coov);

    function plot_coov_density(c)
        lc = log10(c);
        %[xi,f] = ksdensity(log10(coov));
        %fill(f,xi,0.6*ones(1,3));
        hist(lc,100);
        pct = prctile(lc, [.1  99.9]);     
        xt = linspace(pct(1), pct(2), 4);
        set(gca, 'XTick', xt);
        set(gca, 'XTickLabel', round(100*10.^xt)/100);
        
        sf = get_standard_figure_font_sizes;
        axis tight;
        axis square
        set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
        set(gca, 'FontSize', sf.axis_tick_labels);
        xlabel('Coefficient of variation', 'FontSize', sf.axis_labels);
        ylabel('Num. genes', 'FontSize', sf.axis_labels);
        plotSave('figures/quality_filter/gene_coov.png');
        close
    end

    function plot_ac_vs_pnz(q,z)
        [~,density,X,Y]=kde2d([ascolumn(z), ascolumn(q)]);
        f = figure;
        sf = get_standard_figure_font_sizes;
        % dens = log10(density); dens(dens == -Inf) = 0;
        [hc, hc] = contourf(X,Y, sqrt(density),100);
        set(hc, 'LineStyle', 'none')
        
        axis([0 1 0 1])
        axis square;
        colormap(1 - redbluecmap); 
        colorbar
        grid on
        
        hold on
        plot([nz_cutoff, nz_cutoff], [corr_cutoff, 1],  '-k', 'LineWidth', 3); 
        plot([0, nz_cutoff], [corr_cutoff, corr_cutoff], '-k', 'LineWidth', 3); 
        
        set(gca, 'XTick', 0:0.2:1);
        set(gca, 'YTick', 0:0.1:1);
        set(gca, 'FontSize', sf.axis_tick_labels);
        xlabel('Prop. zero expressed', 'FontSize', sf.axis_labels);
        ylabel('Avg. |PCC| with other samples', 'FontSize', sf.axis_labels);
        plotSave('figures/quality_filter/avg_corr_vs_prop_zero_expressed.png');
        close;
    end


    function plot_mr_vs_md(d,mr)
        [~,density,X,Y]=kde2d([log10(ascolumn(d)), ascolumn(mr)]);
        f = figure;
        sf = get_standard_figure_font_sizes;
        [hc, hc] = contourf(X,Y,density,100);
        set(hc, 'LineStyle', 'none')
        
        xmin = 3;
        xmax = 8.3;
        axis([xmin xmax 0 1])
        axis square;
        colormap(1 - redbluecmap); 
        colorbar
        grid on
        
        hold on
        plot(log10([md_cutoff, md_cutoff]), [mr_cutoff, 1], '-k', 'LineWidth', 3); 
        plot(log10([md_cutoff, 10^xmax]), [mr_cutoff, mr_cutoff], '-k', 'LineWidth', 3); 
        
        set(gca, 'XTick', xmin:xmax);
        set(gca, 'YTick', 0:0.1:1);
        set(gca, 'XTickLabel', strcat('10^', cellstr(num2str([3:9]'))'));
        set(gca, 'FontSize', sf.axis_tick_labels);
        xlabel(gca, 'Mapped depth', 'FontSize', sf.axis_labels);
        ylabel(gca, 'Mapped proportion', 'FontSize', sf.axis_labels);
        plotSave('figures/quality_filter/mapped_rate_vs_mapped_depth.png');
        close;
    end


end


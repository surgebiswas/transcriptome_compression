function  lsm = NCBI_SRA_Athalianal_COV_density_plot( Y )


    sig = std(Y,0,2);
    mu = mean(Y,2);
    lsm = log10(sig./mu);
    xt = -0.5:0.5:2;
    xtl = 10.^xt;
    
    
    [ff,xi] = ksdensity(lsm);
    figure;
    jbfill(xi, ff, zeros(1,length(ff)), 'k', 'k', true, 0.5);
    axis square
    set(gca, 'FontSize', 12);
    xlabel('Expression COV (s.d./mean)', 'FontSize', 14);
    ylabel('Density', 'FontSize', 14);
    set(gca, 'XTick', xt);
    set(gca, 'XTickLabel', xtl);
    axis tight;
    axis square;
    plotSave('figures/COV_density.png');
    close


end


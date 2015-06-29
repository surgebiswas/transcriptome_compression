function NCBI_SRA_Athaliana_plot_clock_genes( lY, tids )

    cg = {'LHY', 'CCA1', 'PRR9', 'TOC1', 'ELF3'};
    atg = {'AT1G01060', 'AT2G46830', 'AT2G46790', 'AT5G61380', 'AT2G25930'};
    AT = [];
    for i = 1 : length(atg)
        AT = [AT, lY(:,strcmpi(tids, atg{i}))];
    end
    
    
    pairs = [{'LHY', 'CCA1'}; {'CCA1', 'TOC1'}; {'ELF3', 'TOC1'}];
    for i = 1 : size(pairs,1)
        i1 = strcmpi(pairs{i,1}, cg);
        i2 = strcmpi(pairs{i,2}, cg);
        
        figure;
        plot(AT(:,i1), AT(:,i2), '.k');
        axis square;
        xlabel(pairs{i,1}, 'FontSize', 14);
        ylabel(pairs{i,2}, 'FontSize', 14);
        set(gca, 'FontSize', 12)
        plotSave(sprintf('figures/clock_gene_examples/%s_%s.png', pairs{i,2}, pairs{i,1}));
        iminvert(sprintf('figures/clock_gene_examples/%s_%s.png', pairs{i,2}, pairs{i,1}));
        close
    end

end


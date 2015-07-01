function NCBI_SRA_Mmusculus_preprocess( s, qt, mainDataFile )

    sf = get_standard_figure_font_sizes;
    s.estNumReads = [];


    Y = s.tpm; 
    sids = s.ids;
    tids = s.transcript_id;
    depth = s.depth;
    mapped_ratio = s.mapped_ratio;
    
    fprintf('Original: Num. samples = %0.0f, Num. features = %0.0f\n', length(sids), length(tids));

    qt = qt(sids,:); % Re order rows of metadata table. 
    
    RCTHRESH = 4e6;
    MRTHRESH = 0.70;
    
    % Basic quality check plots. 
    % Mapped percentage vs sample depth.
    if true
        semilogx((s.depth), 100*(s.mapped_ratio), '.k');
        set(gca, 'FontSize', sf.axis_tick_labels);
        xlabel('Mapped depth', 'FontSize', sf.axis_labels);
        ylabel('Mapped percentage', 'FontSize', sf.axis_labels);
        axis square
        hold on
        plot([RCTHRESH RCTHRESH], 100*[MRTHRESH 1], '-r', 'LineWidth', 2);
        plot([RCTHRESH 1e10], 100*[MRTHRESH MRTHRESH], '-r', 'LineWidth', 2);

        v = axis;
        v(2) = 2e8;
        v(1) = 100;
        axis(v);
        set(gca, 'XTick', round(logspace(2,8,7)));

        plotSave('figures/preprocess/mapped_pct_vs_depth.png');
        close;

    end

    % Sample thresholding.
    if true
        torm = s.depth<RCTHRESH & s.mapped_ratio<MRTHRESH;
        Y(:, torm) = [];
        sids(torm) = [];
        qt(torm,:) = [];
        depth(torm) = [];
        mapped_ratio(torm) = [];
    end
    
    fprintf('After mapped pct/depth threshold: Num. samples = %0.0f, Num. features = %0.0f\n', length(sids), length(tids));
    
    
    % Collapse isoform table to be a gene table.
    [tids, Y] = collapse_mouse_isoform_table(tids, Y);
    fprintf('After collapsing isoform table: Num. samples = %0.0f, Num. features = %0.0f\n', length(sids), length(tids));
    
    % TPM profile check
    % TPM profile checking
    if true
        
        yp = TPM_profile_check(Y);
        [yp, pct, torm_pct] = TPM_profile_check(Y, 'zeropropflag',mean(mean(yp == 0)) + 2*std(mean(yp == 0)));
        figure;
        semilogx(yp, pct, '-k');
        hold on
        semilogx(yp(:,torm_pct), pct, '-r');
        axis square
        set(gca, 'FontSize', sf.axis_tick_labels);
        xlabel('TPM', 'FontSize', sf.axis_labels);
        ylabel('Percentile', 'FontSize', sf.axis_labels);
        axis tight;
        axis square
        set(gca, 'XTick', logspace(-2, 5, 8));
        plotSave('figures/preprocess/TPM_profile_check.png');
        close;
        
        Y(:, torm_pct) = [];
        sids(torm_pct) = [];
        depth(torm_pct) = [];
        mapped_ratio(torm_pct) = [];
        qt(torm_pct,:) = [];
    end
    
    
    % Remove lowly expressed genes.
    if true
        lg = mean(Y,2) < 1;
        Y(lg,:) = [];
        tids(lg) = [];
    end
    
    
    save(mainDataFile, 'Y', 'tids', 'sids', 'depth', 'mapped_ratio', 'qt'); 
    fprintf('Output saved in %s\n', mainDataFile);

 

end


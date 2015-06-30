function NCBI_SRA_Athaliana_preprocess(mainDataFile, queryTable)

    sf = get_standard_figure_font_sizes;

    % Read the query table.
    if true
        qt = read_ncbi_sra_query_table(queryTable);
        save([queryTable, '_TEMP_DO_NOT_USE.mat'], 'qt')
    else
        load([queryTable, '_TEMP_DO_NOT_USE.mat']);
    end
    
    % Read TPM data
    load NCBI_SRA_Athaliana_full_data_up_to_18May2015.mat;
    Y = s.tpm; 
    sids = s.ids;
    tids = s.transcript_id;

    qt = qt(sids,:); % Re order rows of metadata table. 
    
    RCTHRESH = 4e6;
    MRTHRESH = 0.75;

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

        plotSave('figures/mapped_pct_vs_depth.png');
        close;

    end

    % Sample thresholding.
    if true
        torm = s.depth<RCTHRESH & s.mapped_ratio<MRTHRESH;
        Y(:, torm) = [];
        sids(torm) = [];
        qt(torm,:) = [];
    end

    % Build the gene count dataset.
    Yd = collapse_Athaliana_isoform_table(mat2dataset(Y, 'ObsNames', tids, 'VarNames', sids));

    % Keep only nuclear protein coding
    dk = dataset('file','/Users/sbiswas/Documents/matlab/src/interactome/At_nuclear_protein_coding.txt', 'ReadObsNames', false, 'ReadVarNames', false);

    o = get(Yd, 'ObsNames');
    ok = dk.Var1;

    [~, ia] = intersect(o, ok);
    npc = setdiff(1:length(o), ia);
    Yd_npc = Yd(npc,:); % Non protein coding features.
    Yd = Yd(ia,:);
    

    Y = double(Yd);
    Y_npc = double(Yd_npc);
    tids = get(Yd, 'ObsNames');
    sids = get(Yd, 'VarNames');
    depth = s.depth;
    mapped_ratio = s.mapped_ratio;
    
    % TPM profile checking
    if true
        
        [yp, pct, torm_pct] = TPM_profile_check(Y);
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
        plotSave('figures/TPM_profile_check.png');
        close;
        
        Y(:, torm_pct) = [];
        Y_npc(:,torm_pct) = [];
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
    
    
    
    % Remove samples with abnormally high non-protein coding features.
    if true
        torm_npc = remove_samples_with_non_prot_coding_contamination(Y, Y_npc, sids);
        
        Y(:, torm_npc) = [];
        Y_npc(:,torm_npc) = [];
        sids(torm_npc) = [];
        depth(torm_npc) = [];
        mapped_ratio(torm_npc) = [];
        qt(torm_npc,:) = [];
    end
    
    

    save(mainDataFile, 'Y', 'tids', 'sids', 'depth', 'mapped_ratio', 'qt'); 
    fprintf('Output saved in %s\n', mainDataFile);



end


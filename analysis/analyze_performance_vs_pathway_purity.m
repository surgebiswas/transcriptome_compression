function analyze_performance_vs_pathway_purity(qt, results, pswap, ktrain, varargin)

csummaries = setParam(varargin, 'compute_summaries', true); % compute summary variables, else make plots.

if csummaries
    pcc_p = zeros(length(results), size(results{1}.s,2));
    pcc_g = zeros(length(results), size(results{1}.z,2));
    uvar_p = pcc_p;
    uvar_g = pcc_g;
    sub = qt.Submission(~ktrain);
    for i = 1 : length(results)
        disp(i);
        sa = subadjust(results{i}.s, sub);
        za = subadjust(results{i}.z, sub);
        sha = subadjust(results{i}.s_hat, sub);
        zha = subadjust(results{i}.z_hat, sub);


        pcc_p(i,:) = rsq_and_slope(sa, sha, 'abs', true); % program sign is arbitrary.
        pcc_g(i,:) = rsq_and_slope(za, zha);

        direc = sign(rsq_and_slope(sa, sha));
        uvar_p(i,:) = var(sa - sha.*repmat(direc,size(sha,1),1))./var(sa);
        uvar_g(i,:) = var(za - zha)./var(za);
    end

    save('perf_vs_pathway_purity_summaries.mat', 'pcc_p', 'pcc_g', 'uvar_p', 'uvar_g');
else   
    
    load('perf_vs_pathway_purity_summaries.mat');
    figure;
    mm = median(pcc_p, 2)';
    L = prctile(pcc_p', 15);
    U = prctile(pcc_p', 85);

    jbfill(pswap, U, L, 'r', 'r', true, 0.3);
    hold on
    h(1) = plot(pswap, mm, '-r', 'LineWidth', 2);

    hold on
    mm = median(pcc_g,2, 'omitnan')';
    L = prctile(pcc_g', 15);
    U = prctile(pcc_g', 85);

    jbfill(pswap, U, L, 'g', 'g', true, 0.3);
    hold on
    h(2) = plot(pswap, mm, '-g', 'LineWidth', 2);
    legend(h, 'tr. programs', 'genes', 'Location', 'SouthWest');
    axis([pswap(end), pswap(1), 0.3 1]);
    axis square

    idxset = [1 2 3 4 5 9];
    set(gca, 'XTick', fliplr(pswap(idxset)))
    set(gca, 'YTick', 0:0.1:1);

    sf = get_standard_figure_font_sizes;
    set(gca, 'FontSize', sf.axis_tick_labels);
    xlabel('Percentage of genes swapped', 'FontSize', sf.axis_labels);
    ylabel('PCC', 'FontSize', sf.axis_labels);
    %set(gca, 'XTick', 'log');
    
    plotSave('figures/perf_vs_pathway_purity/pcc_vs_pathway_purity.png');
    close
    
    % --------- %
    figure;
    mm = median(uvar_p, 2)';
    L = prctile(uvar_p', 15);
    U = prctile(uvar_p', 85);

    jbfill(pswap, U, L, 'r', 'r', true, 0.3);
    hold on
    h(1) = plot(pswap, mm, '-r', 'LineWidth', 2);

    hold on
    mm = median(uvar_g,2, 'omitnan')';
    L = prctile(uvar_g', 15);
    U = prctile(uvar_g', 85);

    jbfill(pswap, U, L, 'g', 'g', true, 0.3);
    hold on
    h(2) = plot(pswap, mm, '-g', 'LineWidth', 2);
    legend(h, 'tr. programs', 'genes', 'Location', 'NorthWest');
    axis([pswap(end), pswap(1), 0 2]);
    axis square

    idxset = [1 2 3 4 5 9];
    set(gca, 'XTick', fliplr(pswap(idxset)))
    set(gca, 'YTick', 0:0.2:2);


    sf = get_standard_figure_font_sizes;
    set(gca, 'FontSize', sf.axis_tick_labels);
    xlabel('Percentage of genes swapped', 'FontSize', sf.axis_labels);
    ylabel('Normalized unexp. var.', 'FontSize', sf.axis_labels);
    
    
    plotSave('figures/perf_vs_pathway_purity/uvar_vs_pathway_purity.png');
    close
end
    
end
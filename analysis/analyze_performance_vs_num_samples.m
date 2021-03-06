function analyze_performance_vs_num_samples(qt, results, nsamples, nsubs, ktrain, varargin)

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
    save('perf_vs_num_samples_summaries.mat', 'pcc_p', 'pcc_g', 'uvar_p', 'uvar_g');
    
else
    
    load('perf_vs_num_samples_summaries.mat');
    figure;
    mm = median(pcc_p,2, 'omitnan')';
    L = prctile(pcc_p', 15);
    U = prctile(pcc_p', 85);

    fh(1) = jbfill(nsamples, U, L, 'r', 'r', true, 0.3);
    hold on
    h(1) = plot(nsamples, mm, '-r', 'LineWidth', 2);

    hold on
    mm = median(pcc_g,2, 'omitnan')';
    L = prctile(pcc_g', 15);
    U = prctile(pcc_g', 85);

    fh(2) = jbfill(nsamples, U, L, 'g', 'g', true, 0.3);
    hold on
    h(2) = plot(nsamples, mm, '-g', 'LineWidth', 2);
    legend(h, 'tr. programs', 'genes', 'Location', 'SouthEast');
    axis([nsamples(1), nsamples(end), 0 1]);
    axis square

    idxset = [1,5, 10, length(nsamples)];
    set(gca, 'XTick', nsamples(idxset))
    set(gca, 'YTick', 0:0.1:1);

    labs = {};
    for i = 1 : length(nsamples)
        labs{i} = [num2str(nsamples(i)), ' (', num2str(nsubs(i)), ')'];
    end
    set(gca, 'XTickLabel', labs(idxset));
    %try; rotateXLabels(gca, 45);end

    sf = get_standard_figure_font_sizes;
    set(gca, 'FontSize', sf.axis_tick_labels);
    xlabel('Training set size', 'FontSize', sf.axis_labels);
    ylabel('PCC', 'FontSize', sf.axis_labels);
    set(gca, 'XScale', 'log');

    plotSave('figures/perf_vs_num_samples/perf_vs_num_samples_pcc.png');






    %%%% -----  %%%%

    figure;
    mm = median(uvar_p,2, 'omitnan')';
    L = prctile(uvar_p', 15);
    U = prctile(uvar_p', 85);

    jbfill(nsamples, U, L, 'r', 'r', true, 0.3);
    hold on
    h(1) = plot(nsamples, mm, '-r', 'LineWidth', 2);

    hold on
    mm = median(uvar_g,2, 'omitnan')';
    L = prctile(uvar_g', 15);
    U = prctile(uvar_g', 85);

    jbfill(nsamples, U, L, 'g', 'g', true, 0.3);
    hold on
    h(2) = plot(nsamples, mm, '-g', 'LineWidth', 2);
    legend(h, 'tr. programs', 'genes', 'Location', 'NorthEast');
    v = axis;
    axis([nsamples(1), nsamples(end), 0 v(4)]);
    axis square

    idxset = [1,5, 10, length(nsamples)];
    set(gca, 'XTick', nsamples(idxset))
    set(gca, 'YTick', 0:0.2:v(4));

    labs = {};
    for i = 1 : length(nsamples)
        labs{i} = [num2str(nsamples(i)), ' (', num2str(nsubs(i)), ')'];
    end
    set(gca, 'XTickLabel', labs(idxset));
     %rotateXLabels(gca, 45);

    sf = get_standard_figure_font_sizes;
    set(gca, 'FontSize', sf.axis_tick_labels);
    xlabel('Training set size', 'FontSize', sf.axis_labels);
    ylabel('Normalized unexp. var.', 'FontSize', sf.axis_labels);
    set(gca, 'XScale', 'log');

    plotSave('figures/perf_vs_num_samples/perf_vs_num_samples_uvar.png');

    close all

end





end


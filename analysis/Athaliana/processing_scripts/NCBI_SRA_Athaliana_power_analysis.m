function NCBI_SRA_Athaliana_power_analysis( lY, qt, latest_training_date, tids, varargin )
%

results = setParam(varargin, 'poweranalysisresults', []);


% See how number of samples correlates with number of submissions rarefied
% to.
if false
% Explore how many samples total are contained in different size subsets
% of submissions.

sub = qt.Submission;
subu = unique(sub);
    N = 10 : 5 : length(subu);
    NREP = 20;
    dm = zeros(NREP, length(N));
    for k = 1 : NREP
        fprintf('%0.0f\n', k);
        sp = generate_submission_subsets(subu, N);
        for j = 1 : length(sp)
            dm(k,j) = sum( steq(sub, sp{j}) );
        end
    end
    
    
    boxplot(dm);
    set(gca, 'XTick', 1:length(N))
    set(gca, 'XTickLabel', N);
    set(gca, 'FontSize', 6);
    set(gca, 'YTick', 0:50:2400);
    xlabel('Number of rarefied submissions', 'FontSize', 14);
    ylabel('Number of samples', 'FontSize', 14);
    plotSave('figures/power_analysis/num_samples_vs_num_submissions.png');
end

% Perform the power analysis.
% We will begin by creating a validation set.
%   e.g. Taking all samples after a certain date.
%
% We will then create a schedule of training sets that vary in size. These
% training set submissions will not overlap with the validation set.
%
% For each training set, we will run marker OMP on it until 75% of the
% variance is explained. 
%
% We will then run tratrain using these selected markers, and will finally
% tradict the validation set.
%
% We will store the following information:
% - global R^2 between tradiction and actual.
% - intra-submission R^2 between tradiction and actual.
% - Actual global and submission adjusted predictions for for 20 random genes. 
% - Marker OMP decomposition results (to assess marker convergence over
%   time).
% - Random generator seed (as a replicate ID).
% - Number of samples/libraries at each submission rarefaction.
%
% Analysis will be run on the PW server.
if false
    % Constants
    LAMBDASCHEDULE = logspace(-2,2,40);
    NSUBMIN = 15;
    SCHEDULESIZE = 15;
    PROPVARUNEXPCUTOFF = 0;
    MAXFEATS = 100;
    rng('default'); KEEPPREDSFOR = randsample(size(lY,2), 50); % always take the same  genes.
    
    SEED = 1; % Change this number to perform replicates
    rng(SEED);
    
    outFile = ['NCBI_SRA_Athaliana_tradict_power_analysis_punexp_', num2str(100*PROPVARUNEXPCUTOFF), ...
        '_maxfeats_', num2str(MAXFEATS), '_seed_', num2str(SEED), '.mat'];
    
    
    % Create the validation set.
    date_cutoff = datenum(latest_training_date);  
    dn = datenum(qt.release_date);
    q = dn >= date_cutoff;
    
    lY_validate = lY(q,:);
    lY_train = lY(~q,:);
    qt_validate = qt(q,:);
    qt_train = qt(~q,:);
    
    
    % Schedule of training set sizes
    sub_train = qt_train.Submission;
    subu_train = unique(sub_train);
    nsub_schedule = round(linspace(NSUBMIN,length(subu_train), SCHEDULESIZE));
    
    
    % Begin the power analysis
    results.seed = SEED;
    results.latest_training_date = latest_training_date;
    results.prop_var_unexp_cutoff = PROPVARUNEXPCUTOFF;
    results.max_feats = MAXFEATS;
    results.lambda_schedule = LAMBDASCHEDULE;
    results.nsubmissions = nsub_schedule;
    results.nsamples = zeros(1, length(nsub_schedule));
    results.somps = cell(1, length(nsub_schedule));
    results.tradict_models = cell(1, length(nsub_schedule));
    results.perfstats = cell(1, length(nsub_schedule));
    
    sp = generate_submission_subsets(subu_train,nsub_schedule);
    for j = 1 : length(sp)
        fprintf('Running iteration %0.0f of %0.0f.\n', j, length(sp));
        
        % Partition out the training set
        ts = steq(sub_train, sp{j});
        lyt = lY_train(ts,:);
        
        % Run marker OMP decomposition.
        somp = marker_OMP(standardize(lyt), PROPVARUNEXPCUTOFF, ... 
                'savememory', true, 'storecrosscorr', true, 'maxfeatures', MAXFEATS);
        lyt_markers = lyt(:,somp.S);
        
        % Tratraining
        model = tratrain(lyt, lyt_markers, 'lambda', LAMBDASCHEDULE);
        
        % Tradicting
        lyv_markers = lY_validate(:,somp.S);
        yhat = tradict(lyv_markers, model);
        
        % Evaluate performance.
        perfstats = evaluate_tradiction(lY_validate, yhat, ...
            'submissionids', qt_validate.Submission, 'keeppredforidx', KEEPPREDSFOR);
       
        results.nsamples(j) = sum(ts);
        results.somps{j} = somp;
        results.tradict_models{j} = model;
        results.perfstats{j} = perfstats;
    end
    
    save(outFile, 'results', '-v7.3');
end


% Graphics of the power analysis
if true
    if isempty(results)
        outFile = ['NCBI_SRA_Athaliana_tradict_power_analysis_punexp_', num2str(100*PROPVARUNEXPCUTOFF), ...
            '_maxfeats_', num2str(MAXFEATS), '_seed_', num2str(SEED), '.mat'];

        load(outFile);
    end
        
    glrsq = [];
    sarsq = [];
    for j = 1 : length(results.perfstats);
        glrsq = [glrsq, results.perfstats{j}.global_Rsq'];
        sarsq = [sarsq, results.perfstats{j}.sub_adj_Rsq'];
    end
    
    glrsq(isnan(sum(glrsq,2)),:) = [];
    sarsq(isnan(sum(sarsq,2)),:) = [];
    
    tk = [1:size(glrsq,2)];
    glxt = 1:3:length(tk)*3;
    saxt = glxt + 1;
    gp = prctile(glrsq, 2.5);
    sp = prctile(sarsq, 2.5);
    
    
    %%% Violin plots
    if false
        
        % Plot Global and Submission Rsq else only global R^2.
        if false
            figure;
            hold on;
            for j = 1 : length(tk)
                g = glrsq(:,tk(j)); g(g < gp(j)) = [];
                s = sarsq(:,tk(j)); s(s < sp(j)) = [];


                violin(g, 'x', glxt(j), 'facecolor', 'm');
                [h, l] = violin(s, 'x', saxt(j), 'facecolor', 'c');
            end
            axis([0.5 max(saxt)+0.5 -0.8 1]);
            grid on
            set(gca, 'XTick', mean([glxt;saxt]));
            set(gca, 'XTickLabel', results.nsubmissions(tk));
            set(l,'Location', 'SouthWest');
            plotSave('figures/power_analysis/violin_plots_Rsq_vs_num_submissions.png');
        else
            figure;
            hold on;
            for j = 1 : length(tk)
                g = glrsq(:,tk(j)); g(g < gp(j)) = [];
                [h, l] = violin(g, 'x', glxt(j), 'facecolor', 'm');
            end
            axis([0.5 max(glxt)+0.5 0 1]);
            grid on
            set(gca, 'XTick', glxt);
            set(gca, 'XTickLabel', results.nsamples(tk));
            set(l,'Location', 'SouthWest');
            plotSave('figures/power_analysis/violin_plots_Rsq_vs_num_submissions_global_Rsq_only.png');
            iminvert('figures/power_analysis/violin_plots_Rsq_vs_num_submissions_global_Rsq_only.png');
        end
    end
    
    %%% Projection plots of improvement.
    if true
        figure;
        projectionPlot(sarsq, results.nsamples, [10, 50, 90], 2, {'r', 'b', 'k'});
        plotSave('figures/power_analysis/submission_adjusted_Rsq_vs_SRA_size_percentiles.png');
        iminvert('figures/power_analysis/submission_adjusted_Rsq_vs_SRA_size_percentiles.png');
        title('Submission adjusted performance', 'FontSize', 16);
        close

        figure;
        projectionPlot(glrsq, results.nsamples, [10, 50, 90], 2, {'r', 'b', 'k'});
        plotSave('figures/power_analysis/global_Rsq_vs_SRA_size_percentiles.png');
        iminvert('figures/power_analysis/global_Rsq_vs_SRA_size_percentiles.png');
        title('Global performance', 'FontSize', 16);
        close
    end
    
    % Marginal utility of each additional, number of markers selected.
    if false
        nm = zeros(length(results.somps),1);
        for j = 1 : length(nm)
            nm(j) = length(results.somps{j}.S);
        end
        subplot(1,2,1)
        plot(results.nsubmissions, nm, '-ok', 'LineWidth', 3);
        set(gca, 'FontSize', 12);
        xlabel('Database size (# submissions)', 'FontSize', 14);
        ylabel('Number of markers selected', 'FontSize', 14);
        
        subplot(1,2,2);
        semilogy(diff(100 - 100*[1 results.somps{end}.punexp]), '-k', 'LineWidth', 3)
        set(gca, 'FontSize', 12);
        xlabel('Marker', 'FontSize', 14);
        ylabel('Marginal contribution to % variance explained', 'FontSize', 14);
        axis tight
        
        plotSave('figures/power_analysis/marker_growth_and_final_marginal_utility.png');
        close
    end
    
    % Plots of example genes, prediction vs. actual.
    if false
        NGENES = 9;
        NROWS = 3;
        NCOLS = 3;
        
        figure;
        for j = 1 : NGENES
            subplot(NROWS, NCOLS,j);
            pred_v_actual_plot(results.perfstats{end}.kept_yhat(:,j), ...
                results.perfstats{end}.kept_ytrue(:,j), ...
                tids{results.perfstats{end}.kept_pred_idx(j)});
        end
        
        plotSave('figures/power_analysis/example_genes_full_training_set_predicting_validation.png');
        iminvert('figures/power_analysis/example_genes_full_training_set_predicting_validation.png');
        close
        
        figure
        hist(results.perfstats{end}.global_Rsq, 100);
        set(gca, 'FontSize', 14);
        axis square;
        vv = axis;
        axis([-0.2 1.01 vv(3) vv(4)]);
        xlabel('R^2 (prediction vs. actual)', 'FontSize', 16);
        ylabel('Number of genes', 'FontSize', 16);
        plotSave('figures/power_analysis/performance_histogram_of_gene_Rsq.png');
        iminvert('figures/power_analysis/performance_histogram_of_gene_Rsq.png');
        close
    end
    
   
end

    function pred_v_actual_plot(x,y,gene_id)
        plot(x,y, '.k');
        axis tight;
        buffer_axis;
        v = axis;
        bc = min([v(1), v(3)]);
        tc = max([v(2), v(4)]);
        
        plot([bc tc], [bc tc], '-r', 'LineWidth', 2);
        hold on
        plot(x,y, '.k'); % put points in front of the line
        axis([bc tc bc tc]);
        title(gene_id);
        axis square;
    end

    function projectionPlot(rsq, nsubmissions, p, burnin, cm)
        rsq(:,1:burnin) = [];
        nsubmissions(1:burnin) = [];
        
        
        pct = prctile(rsq,p)';
        
        for i = 1 : size(pct,2)
            x = [ones(length(nsubmissions),1), log10(nsubmissions')];
            b = pinv(x)*pct(:,i);

            semilogx(nsubmissions, pct(:,i), 'o', 'Color', cm{i}, 'MarkerSize', 12, 'LineWidth', 3)
            hold on
            nsub_extrap = 4*nsubmissions(end);
            xquery = [x; [1 log10(nsub_extrap)]];
            h(i) = semilogx([nsubmissions, nsub_extrap], xquery*b, '-', 'Color', cm{i}, 'LineWidth', 3);
        end
        
        xt = [[0.25 0.5 0.75 1]*max(nsubmissions), nsub_extrap];
        %xt = [logspace(log10(min(nsubmissions)), log10(max(nsubmissions)), 4), nsub_extrap];
        set(gca, 'XTick', round(xt));
        set(gca, 'XTickLabel', round(xt));
        axis square;
        axis tight;
        set(gca, 'FontSize', 12);
        
        lv = [];
        for i = 1 : length(p)
            lv = [lv, {[num2str(p(i)), 'th %-ile']}];
        end
        l = legend(h, lv, 'Location', 'SouthEast');
        set(l, 'Box', 'off');
        ylabel('R^2 (prediction vs. actual)', 'FontSize', 14);
        xlabel('Database size (# submissions)', 'FontSize', 14);
        
    end



    function sp = generate_submission_subsets(subu,n)
        sp = cell(1, length(n));
        for i = 1 : length(n)
            ind = randsample(length(subu), n(i));
            sp{i} = subu(ind);
        end
        
    end

end


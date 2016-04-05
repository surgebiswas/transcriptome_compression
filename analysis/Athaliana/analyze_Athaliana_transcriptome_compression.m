clear;
rng('default');

% Initialize environment
completiontext = false;
[~,hn] = unix('hostname'); % Get hostname.
if ~isempty(strfind(hn, '.kd.unc.edu'))
    % We are on killdevil
    homedir = '/proj/dangl_lab/sbiswas';
    completiontext = true;
elseif strcmpi(strtrim(hn), 'ygritte')
    homedir = '/home/sbiswas';
else 
    % We are working locally
    homedir = '/Users/sbiswas';
    completiontext = false;
end
repo = 'transcriptome_compression/';
datadir = [homedir, '/GitHub/data/', repo, 'Athaliana/'];
path(genpath([homedir, '/GitHub/', repo]), path);
cd(datadir);
tic;

%mainDataFile = 'NCBI_SRA_Athaliana_full_data_up_to_18May2015_processed.mat';
mainDataFile = 'NCBI_SRA_Athaliana_full_data_up_to_06Sept2015_quality_filtered.mat';
queryTable = 'Athaliana_query_table_06Sept2015.csv';

% Update raw (unquality filtered) collection of transcriptomes.
if false
    load('NCBI_SRA_Athaliana_full_data_up_to_18May2015.mat');
    [~,sidx] = sort(s.transcript_id);
    s.tpm = s.tpm(sidx,:);
    s.transcript_id = s.transcript_id(sidx);
    s1 = s;
    
    load('Athaliana_download_06Sept2015_assembled_srafish_output.mat');
    [~,sidx] = sort(s.transcript_id);
    s.tpm = s.tpm(sidx,:);
    s.transcript_id = s.transcript_id(sidx);
    
    % Ensure transcript ID's are in the same orders. TPM tables are sorted
    % accordingly as well.
    for i = 1 : length(s.transcript_id)
        assert(strcmpi(s.transcript_id{i}, s1.transcript_id{i}));
    end
    
    
    sjoin.ids = [s1.ids, s.ids];
    sjoin.mapped_ratio = [s1.mapped_ratio, s.mapped_ratio];
    sjoin.tpm = [s1.tpm, s.tpm];
    sjoin.transcript_id = s.transcript_id; % Known to be the same as s1.
    sjoin.depth = [s1.depth, s.depth];
    
    save('NCBI_SRA_Athaliana_full_data_up_to_06Sept2015.mat', 'sjoin');
end

% Output list of successfully downloaded samples.
if false; 
    fts = regexprep(filetimestamp, '_.+', '');
    
    output_list_of_successful_downloads('NCBI_SRA_Athaliana_full_data_up_to_18May2015.mat', ['NCBI_SRA_Athaliana_successfully_downloaded_', fts, '.txt']);
end

% Quality filter. Sample and expression filtering.
% Isoform collapsing to genes
% Keeping of nuclear protein coding genes.
% Check TPM profiles
if false; 
    qfparams.MRTHRESH = 0.75;
    qfparams.RCTHRESH = 4e6;
    qfparams.CORRCUTOFF = 0.45;
    qfparams.NZCUTOFF = 0.2; 
    qfparams.TPMCUTOFF = 1;
    
    if true
        load('NCBI_SRA_Athaliana_full_data_up_to_06Sept2015.mat');
        [Y, sids, tids] = quality_filter(sjoin, qfparams, 'Athaliana');
    end
    
    if false
        plot_quality_filter_report('Athaliana_quality_filter_summary_data.mat', qfparams, 'Athaliana');
    end

    
    % Prepare the query table
    qt_full = read_ncbi_sra_query_table(queryTable);
    
    % Some entries from the SRA have been removed before we updated
    % the query table. Remove these samples from analysis.
    k = steq(sids, get(qt_full, 'ObsNames'));
    Y = Y(:,k);
    sids = sids(k);
    
    qt = qt_full(sids,:);
    save(mainDataFile, 'Y', 'tids', 'sids', 'qt');

    return; 
end

load(mainDataFile);
lY = log10(Y' + 0.1);

% Number of HQ transcriptomes in the SRA as a function of time.
if false
    plot_ncbi_sra_growth(qt);
    plotSave('figures/growth_over_time/NCBI_SRA_Athaliana_growth_v_time.png');
    close
end


% Plots of Clock genes.
if false; 
    NCBI_SRA_Athaliana_plot_clock_genes(lY, tids); 
end


% Coefficient of variation density plot.
if false; 
    logcov = NCBI_SRA_Athalianal_COV_density_plot(Y); 
end

% Perc. variation explained vs eigengene.
if false; 
    [coef, pexp] = pexp_vs_components(lY, 'Athaliana'); 
end

% PCA plots for main figures
if false
    load('NCBI_SRA_Athaliana_PCA_pexp_vs_eigengene_params.mat');
    NCBI_SRA_Athaliana_plot_PCA( lY, coef, qt, pexp )
end

% PCA plots by submission, for supplemental figures
if false;
    load('NCBI_SRA_Athaliana_PCA_pexp_vs_eigengene_params.mat');
    pca_by_submission(lY, qt, coef, pexp);
    plotSave('figures/pca/NCBI_SRA_Athaliana_PCA_by_submission.png');
    close;
end

% PC Stability plots for supplemental
if false
    load('NCBI_SRA_Athaliana_PCA_pexp_vs_eigengene_params.mat');
    PC_stability_plots(coef,pexp,lY,qt, 'makersqplots', true);
end


% Convergence of data clusters along the principal components.
if false; 
    save('PCA_pexp_vs_eigengene_params.mat', 'coef', 'pexp');
    NCBI_SRA_PC_stability(coef, pexp, lY, qt);
end

% Marker OMP decomposition on FULL data
% Run on $pw
if false 
    % Last run Mar 2, 2016.
    rng('default');
    punexp = 0;
    maxfeats = 500;
    somp = marker_OMP(standardize(lY), punexp, 'savememory', false, 'maxfeatures', maxfeats, 'subresidual', 0.1);
    save(sprintf('NCBI_SRA_Athaliana_marker_OMP_decomposition_punexp_%0.2f_maxfeats_%0.0f.mat', punexp, maxfeats), 'somp');
end

% Compress the FULL data.
% Used for figure 2.
if false
    NMARKERS = 100;
    load('NCBI_SRA_Athaliana_marker_OMP_decomposition_punexp_0.00_maxfeats_500.mat');
    model = tratrain(lY, lY(:, somp.S(1:NMARKERS)));
    model.reconstruction = lY(:, somp.S(1:NMARKERS))*model.b + repmat(model.b0,size(lY,1),1);
    save(['NCBI_SRA_Athaliana_compression_and_reconstruction_nmarkers_', num2str(NMARKERS), '.mat'], 'model')
end


% Preliminary investigation into locally weighted regression
if false
    load NCBI_SRA_Athaliana_marker_OMP_decomposition_punexp_0.00_maxfeats_500.mat;
    nmarkers = 1;
    markers = somp.S(1:nmarkers);
    method = 'ridgefit';
    
    % Loss function
    trimmean_percent = 10;
    lf = @(ytest,yhat) trimmean(abs((ytest(:) - yhat(:))./ytest(:)), ...
        trimmean_percent, 'round', 1); % Robust relative error.
    
    % Cross validation indices.
    % Group by submission.
    nfolds = 10; %length(unique(qt.Submission));
    cvi = kfoldcrossvalindbygroup(nfolds, qt.Submission);
    
    % Representative residual indices.
    rng('default');
    rep_residual = randsample(size(lY,2), round(0.05*size(lY,2)));
    
    if strcmpi(method, 'kernelfit')
        % Kernel fit train parameters
        params_tune{1}{1} = @gaussian_kernel_diagonal; % Kernel
        %params_tune{2} = ([0.1:0.02:0.5]')*sqrt(nmarkers); % loc avg; scale by the dimension of the query space
        params_tune{2} = 2.3*sqrt(nmarkers); % loc reg;
        params_tune{3} = true; % train standardized?
        params_tune{4} = true; % predict original?
        params_tune{5}{1} = 'local_regression'; % ['local_average' | 'local_regression'];

        trainfun = @kernelfit_train;
        predfun = @kernelfit_predict;
    end
    if strcmpi(method, 'ridgefit')
        params_tune{1} = [0, logspace(-6,-2,19)]';
        params_tune{2} = true;
        
        % Precompute x'*x and x'*y
        params_tune{3}{1} = ridgefit_precompute;
        params_tune{3}{1}.itr = 1;
        params_tune{3}{1}.precompute = cell(nfolds,2);
        for i = 1 : nfolds
            x = lY(cvi~=i,markers); x = [ones(size(x,1),1),x];
            y = lY(cvi~=i,rep_residual);
            
            params_tune{3}{1}.precompute{i,1} =x'*x;
            params_tune{3}{1}.precompute{i,2} =x'*y;
        end
        
        trainfun = @ridgefit_train;
        predfun = @ridgefit_predict;
    end
    
    
    %profile on;
    cvs = cvalidate_tune(lY(:,markers), lY(:,rep_residual), trainfun, ...
        predfun, params_tune, 'crossvalind', cvi, 'loss_fun', lf);
    %profile off;
    
    %profile viewer
end


% Preliminary investigation into compression ensembles.
if false
    rng('default');
    [ytrain, ytest, trainind] = partition_data(lY, qt, 0.1);
    
    if false
    compression_ensemble_train(ytrain, 'saveFile', 'NCBI_SRA_Athaliana_compression_ensemble_test.mat');
    end
    
    if false
        % Let's compare this to a 500 marker global decomposition.
        rng('default');
        gl_stats = marker_OMP(standardize(ytrain), 0, 'maxfeatures', 500, 'subresidual', 0.01);
        save('NCBI_SRA_Athaliana_compression_ensemble_test.mat', 'gl_stats', '-append');
    end
end


% Heatmap of original vs. reconstruction of FULL data.
% Part of Figure 2.
if false
    NMARKERS = 100;
    load('NCBI_SRA_Athaliana_marker_OMP_decomposition_punexp_0.00_maxfeats_500.mat');
    load(['NCBI_SRA_Athaliana_compression_and_reconstruction_nmarkers_', num2str(NMARKERS), '.mat'])
    
    model.reconstruction = tradict( lY(:,somp.S(1:NMARKERS)), model);
    
    heatmap_raw_vs_reconstructed(Y',somp, model, 'Athaliana', false, NMARKERS);
    
    f = pred_v_actual_density_plot(lY, model.reconstruction, 'subsample_genes', 5000);
    hcb=colorbar;
    ca = caxis;
    set(hcb,'YTick',round(100*linspace(min(ca), max(ca), 4))/100);
    sf = get_standard_figure_font_sizes;
    set(gca, 'FontSize', sf.axis_tick_labels);
    plotSave('figures/heatmap_original_vs_reconstruction/insample_density_plot.png');
    close
end

% Gene set analysis.
if true
    load('~/GitHub/transcriptome_compression/analysis/gene_ontology/Athaliana_representative_gene_set_02-Apr-2016.mat');
%     load('NCBI_SRA_Athaliana_cluster_idx_for_raw_vs_reconstructed_heatmap.mat');
%     [sY, train_mu, train_sig] = standardize(lY);
%     [ys, ngenes, engenes, pexp, gscoef] = collapse_to_gene_sets(sY, tids, sets);
%     save('NCBI_SRA_Athaliana_gene_set_PC_coefs.mat', 'ngenes', 'engenes', 'pexp', 'gscoef', 'train_mu', 'train_sig');
    
    if false
        rng('default')
        [sY, stats.train_mu, stats.train_sig] = standardize(lY);
        stats = geneset_cluster( sY, tids, sets, 'stats', stats );
        stats = geneset_encode(sY, 100, stats);

        save('NCBI_SRA_Athaliana_geneset_encoding.mat', 'stats');
    end
    
    
    % Ridge fit training for predicting gene sets.
    if false
        load('NCBI_SRA_Athaliana_geneset_encoding.mat');
        markers = stats.S;
        nfolds = length(unique(qt.Submission));
        cvi = kfoldcrossvalindbygroup(nfolds, qt.Submission);
        target = stats.geneset.sy_sets;
        
        
        % Loss function
        trimmean_percent = 10;
        %lf = @(ytest,yhat) sqrt( var(ytest(:)-yhat(:),1) );
        lf = @(ytest,yhat) trimmean(abs((ytest(:) - yhat(:))./ytest(:)), ...
            trimmean_percent, 'round', 1); % Robust relative error.
        
        
        params_tune{1} = [0, logspace(-6,-2,19)]';
        params_tune{2} = true;
        
        % Precompute x'*x and x'*y
        params_tune{3}{1} = ridgefit_precompute;
        params_tune{3}{1}.itr = 1;
        params_tune{3}{1}.precompute = cell(nfolds,2);
        for i = 1 : nfolds
            x = lY(cvi~=i,markers); x = [ones(size(x,1),1),x];
            y = target(cvi~=i,:);
            
            params_tune{3}{1}.precompute{i,1} =x'*x;
            params_tune{3}{1}.precompute{i,2} =x'*y;
        end
        
        trainfun = @ridgefit_train;
        predfun = @ridgefit_predict;
        
        
        % Cross validate to select L2 tuning parameter.
        cvs = cvalidate_tune(lY(:,markers), target, trainfun, ...
        predfun, params_tune, 'crossvalind', cvi, 'loss_fun', lf);
        
    
        save('NCBI_SRA_Athaliana_geneset_encoding.mat', 'stats', 'cvs');
    end
    
    if true
        load('NCBI_SRA_Athaliana_geneset_encoding.mat');
        target = stats.geneset.sy_sets;
        markers = stats.S;
        trainfun = @ridgefit_train;
        predfun = @ridgefit_predict;
        ps = evaluate_prospective_performance_2(target,lY(:,markers), qt, predfun, trainfun, cvs.param_best_perf);
        plotSave('figures/prospective_performance/losocv_genesets_density.png');
        close
        
        lYs = standardize(ps.lY);
        lYsh = standardize(ps.lY_hat);
        [ri,ci] = hclust(lYs);
        
        imagesc(lYs(ri,ci),[-3 3]); colormap(prgn)
        axis off
        plotSave('figures/prospective_performance/losocv_genesets_heatmap_actual.png');
        close
        
        imagesc(lYsh(ri,ci),[-3 3]); colormap(prgn)
        axis off
        plotSave('figures/prospective_performance/losocv_genesets_heatmap_predicted.png');
        close

    end
    
    
    
    
end

%%% PROSPECTIVE PERFORMANCE
% Train on 90% of data (cutoff determined by date)
% Test on remaining 10%. 
if false
    if false
        evaluate_prospective_performance(lY,qt, 'NCBI_SRA_Athaliana_prospective_performance.mat');
    else
        load('NCBI_SRA_Athaliana_prospective_performance.mat');
        [ perfstats ] = evaluate_tradiction( ytest, yhat, 'submissionids', qt.Submission(~trainind));
        
        set(0, 'currentfigure', perfstats.figs.fh_density_global);
        plotSave('figures/prospective_performance/global_density_training_data_100pct.png')
        set(0, 'currentfigure', perfstats.figs.fh_density_subadj);
        plotSave('figures/prospective_performance/intrasubmission_density_training_data_100pct.png')
        
        set(0, 'currentfigure', perfstats.figs.fh_hist_global);
        plotSave('figures/prospective_performance/global_histogram_training_data_100pct.png')
        set(0, 'currentfigure', perfstats.figs.fh_hist_subadj);
        plotSave('figures/prospective_performance/intrasubmission_histogram_training_data_100pct.png')
        
        set(0, 'currentfigure', perfstats.figs.cbar);
        plotSave('figures/prospective_performance/colorbar.png');
        close all
        
        
    end
end

% Test of residual subsampling. 
% Compare to result of compressing full data.
if false
    punexp = 0;
    maxfeats = 100;
    subsampleto = [0.01 0.05 0.1];
    
    for i = 1 : length(subsampleto)
        somp{i} = marker_OMP(standardize(lY), punexp, 'savememory', true, 'maxfeatures', maxfeats, 'subresidual', subsampleto(i));
        save(sprintf('NCBI_SRA_Athaliana_marker_OMP_decomposition_residual_subsample_punexp_%0.2f_maxfeats_%0.0f.mat', punexp, maxfeats), 'somp', 'subsampleto');
    end
end


% PURGATORY
if false

    % Tradict new data from old
    if false
        load('NCBI_SRA_Athaliana_marker_OMP_decomposition.mat');

        NCBI_SRA_Athaliana_tradict_new_from_old( lY, qt, somp, '05-Mar-2015' );

    end

    % Tradiction power analysis
    if false
        if true
            % perform the analysis and make graphics
            PROPVARUNEXPCUTOFF = 0;
            MAXFEATS = 100;
            SEED = 1;
            matfile = ['NCBI_SRA_Athaliana_tradict_power_analysis_punexp_', num2str(100*PROPVARUNEXPCUTOFF), ...
                    '_maxfeats_', num2str(MAXFEATS), '_seed_', num2str(SEED), '.mat'];

            load(matfile) % results in memory.
            NCBI_SRA_Athaliana_power_analysis(lY,qt, '05-Mar-2015', tids, 'poweranalysisresults', results);
        else
            NCBI_SRA_Athaliana_power_analysis(lY,qt, '05-Mar-2015', tids);
        end
    end

    % Context specific performance failed.
    if false
        NCBI_SRA_Athaliana_context_specific_perfomance(lY,qt, '05-Mar-2015');
    end



    % Supervised model. Tradiction. OLD CODE (see above).
    if false

        if false
            load('PCA_pexp_vs_eigengene_params.mat');
            NMARKERS = 100;
            I = eye(100);
            Rs = corr( coef(:,1:size(I,2))', I);

            [~,mind] = max(Rs);

            cvind = crossvalind('Kfold',size(lY,1),10);

            lY_train = lY(cvind ~= 1,:);
            x_train = lY_train(:,mind);

            model = tratrain(lY_train,x_train,'lambda', 55);

            lY_validate = lY(cvind == 1, :);
            x_validate = lY_validate(:,mind);

            yhat = tradict(x_validate, model);

            Rsq = zeros(1, size(yhat,2));
            slopes = zeros(1, size(yhat,2));
            for i = 1 : length(Rsq)
                Rsq(i) = corr(yhat(:,i), lY_validate(:,i));
                c = pca([lY_validate(:,i), yhat(:,i)]);
                slopes(i) = c(2,1)/c(1,1);
            end

            save('tradict_results.mat');
        else
            load('tradict_results.mat');

            subplot(2,1,1)
            hist(Rsq,100);
            axis square
            ylabel('Frequency');
            xlabel('R^2 (prediction vs. actual)');

            subplot(2,1,2)
            hist(slopes,100);
            axis square
            ylabel('Frequency');
            xlabel('Slope (prediction vs. actual)');
            plotSave('figures/test_set_overall_prediction_quality.png');
            iminvert('figures/test_set_overall_prediction_quality.png');
            close

            rind = randsample(length(Rsq), 8, false);
            for i = 1 : length(rind)
                subplot(2,4,i)

                plot(lY_validate(:,rind(i)), yhat(:,rind(i)), '.k')
                title(tids{rind(i)});
                axis square;
                axis tight;
                buffer_axis;

                v = axis;
                v(1) = v(3);
                v(2) = v(4);
                axis(v);

                hold on
                plot([v(1) v(2)], [v(3) v(4)], '-r', 'LineWidth', 2);
            end
            plotSave('figures/test_set_prediction_examples.png');
            iminvert('figures/test_set_prediction_examples.png');
            close

        end

    end

    % True test. Comparison to Col-0 temperature perturbation.
    if false
        old = cd('/Users/sbiswas/Documents/matlab/src/tradict/NCBI_SRA/Col0_temperature_test');
        dk = dataset('file','/Users/sbiswas/Documents/matlab/src/interactome/At_nuclear_protein_coding.txt', 'ReadObsNames', false, 'ReadVarNames', false);
        dko = dk.Var1;

        s22 = read_sailfish_output('Col0_22C_ZT0_quant_bias_corrected.sf');
        s27 = read_sailfish_output('Col0_27C_ZT0_quant_bias_corrected.sf');
        s22c = collapse_isoform_table(s22.table);
        s27c = collapse_isoform_table(s27.table);

        d22 = dataset('file', 'Col0_22C_ZT0.ct', 'ReadObsNames', true, 'ReadVarNames', false);
        d27 = dataset('file', 'Col0_27C_ZT0.ct', 'ReadObsNames', true, 'ReadVarNames', false);

        % Keep common features (should only be nuclear protein coding).
        [~,ia] = intersect(get(s22c, 'ObsNames'), tids);
        s22c = s22c(ia,:);
        s27c = s27c(ia,:);

        [~,ia] = intersect(get(d22, 'ObsNames'), tids);
        d22 = d22(ia,:);
        d27 = d27(ia,:);


        st = log10([s22c.TPM, s27c.TPM]+0.1); 
        dt = log10(1e6*relabund([d22.Var2, d27.Var2]'+1)');


        slfc = st(:,2) - st(:,1);
        dlfc = dt(:,2) - dt(:,1);    


        % Evaluate prediction accuracy.
        load('tradict_results.mat');


        yhat_col = tradict(st(mind,:)', model)';
        slfc_pred = yhat_col(:,2) - yhat_col(:,1);

        figure
        plottok = prod([d22.Var2, d27.Var2],2) > 0;
        plotmatrix([dt(plottok,1), st(plottok,1), yhat_col(plottok,1)], '.k')
        plotSave('figures/plotmatrix_standard_sailfish_sailfishpred_22C.png');
        iminvert('figures/plotmatrix_standard_sailfish_sailfishpred_22C.png');

        plotmatrix([dt(plottok,2), st(plottok,2), yhat_col(plottok,2)], '.k')
        plotSave('figures/plotmatrix_standard_sailfish_sailfishpred_27C.png');
        iminvert('figures/plotmatrix_standard_sailfish_sailfishpred_27C.png');

        close all


        figure;
        plot(dt(plottok,1), dt(plottok,2), '.k');
        axis square
        xlabel(['Col-0 22C ZT0', 10, 'Standard log_{10}(CPM)']);
        ylabel(['Col-0 27C ZT0', 10, 'Standard log_{10}(CPM)']);
        plotSave('figures/Col0_ZT0_27C_vs_22C_standard_pipeline.png');
        iminvert('figures/Col0_ZT0_27C_vs_22C_standard_pipeline.png');
        close


        figure;
        [ff,xi] = ksdensity(dt(:,2) - dt(:,1));
        jbfill(xi, ff, zeros(1,length(ff)), [0.3 0.3 0.3], [0.3 0.3 0.3], true, 0.6);
        axis square
        xlabel('log_2(27C/22C) Col-0 ZT0');
        ylabel('Density');
        plotSave('figures/Col0_ZT0_27C_vs_22C_standard_pipeline_log2FC.png');
        iminvert('figures/Col0_ZT0_27C_vs_22C_standard_pipeline_log2FC.png');
        close


        figure
        plot(dlfc(plottok), slfc(plottok), '.k');
        text(-1.5, 1.5, sprintf('R^2 = %0.2f', corr(dlfc(plottok), slfc(plottok))));
        axis([-2 2 -2 2]);
        axis square
        xlabel('log_{10}(27C/22C) Col-0 ZT0 Standard Pipeline');
        ylabel('log_{10}(27C/22C) Col-0 ZT0 Sailfish');
        set(gca, 'FontSize', 14)
        plotSave('figures/logFC_standard_vs_sailfish.png');
        iminvert('figures/logFC_standard_vs_sailfish.png');
        close

        figure
        plot(dlfc(plottok), slfc_pred(plottok), '.k');
        text(-1.5, 1.5, sprintf('R^2 = %0.2f', corr(dlfc(plottok), slfc_pred(plottok))));
        axis square
        axis([-2 2 -2 2]);
        xlabel('log_{10}(27C/22C) Col-0 ZT0 Standard Pipeline');
        ylabel('log_{10}(27C/22C) Col-0 ZT0 Tradict');
        set(gca, 'FontSize', 14)
        plotSave('figures/logFC_standard_vs_tradict.png');
        iminvert('figures/logFC_standard_vs_tradict.png');
        close



        cd(old)
    end


    % Toy example of tradict
    if false
        n = 67;
        p = 100;
        x = zeros(n, p);
        kind = [25 80];

        x(:, kind(1)) = randn(n,1);
        x(:, kind(2)) = randn(n,1);
        [~,sidx] = sort(x(:,kind(1)));

        b = [[0.7*randn(1, p-72)-1; 0.7*randn(1,p-72) + 2], [0.7*randn(1, p-30) + 2; 0.7*randn(1,p-30) - 1]]  ;
        xfull = x;
        xfull(:, setdiff(1:p,kind)) = x(:,kind)*b;

        figure;
        imagesc(x(sidx,kind), [-2 2]); colormap(prgn);
        axis image
        set(gca, 'TickLength', [0 0]);
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        plotSave('figures/toy_example/markers_train.png');
        iminvert('figures/toy_example/markers_train.png');
        close

        figure;
        imagesc(standardize(xfull(sidx,setdiff(1:p,kind))), [-2 2]); colormap(prgn)
        axis image
        set(gca, 'TickLength', [0 0]);
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        plotSave('figures/toy_example/full_train.png');
        iminvert('figures/toy_example/full_train.png');
        close

        figure;
        imagesc(b); colormap(redbluecmap)
        axis image
        set(gca, 'TickLength', [0 0]);
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        plotSave('figures/toy_example/coefs.png');
        iminvert('figures/toy_example/coefs.png');
        close

        nnew = 15;
        xnew = randn(nnew,2);
        ynew = xnew*b;

        [~,sidx2] = sort(xnew(:,1));

        imagesc(standardize(xnew(sidx2,:)), [-2 2]); colormap(prgn)
        axis image
        set(gca, 'TickLength', [0 0]);
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        plotSave('figures/toy_example/markers_test.png');
        iminvert('figures/toy_example/markers_test.png');
        close

        imagesc(standardize(ynew(sidx2,:)), [-2 2]); colormap(prgn)
        axis image
        set(gca, 'TickLength', [0 0]);
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        plotSave('figures/toy_example/full_test.png');
        iminvert('figures/toy_example/full_test.png');
        close








    end
end



t = toc;
if completiontext
    msg = sprintf('Job took %0.2f seconds.\n... you''re awesome.', t);
    send_text_message('919-757-1609', 'AT&T', 'Job Finished', msg)
end
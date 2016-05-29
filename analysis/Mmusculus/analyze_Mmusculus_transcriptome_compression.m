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
end
repo = 'transcriptome_compression/';
datadir = [homedir, '/GitHub/data/', repo, 'Mmusculus/'];
path(genpath([homedir, '/GitHub/', repo]), path);
cd(datadir);
tic;

% General variables
qtfile = 'Mmusculus_query_table_19Sept2015.csv';
mainDataFile = 'NCBI_SRA_Mmusculus_full_data_up_to_19Sept2015_quality_filtered.mat';

% Construction full collection of raw transcriptomes.
if false
    load('assembled_srafish_output/NCBI_SRA_Mmusculus_download_04June2015_prelim.mat');
    s1 = s;
    
    load('assembled_srafish_output/Mmusculus_download_11July2015_assembled_srafish_output.mat');
    s2 = s;
    
    load('assembled_srafish_output/Mmusculus_download_19Sept2015_assembled_srafish_output.mat');
    s3 = s;
    
    s = update_raw_transcriptome_collection(s1,s2);
    s = update_raw_transcriptome_collection(s, s3);
    
    save('NCBI_SRA_Mmusculus_full_data_up_to_19Sept2015.mat', 's', '-v7.3');
    return;
end


% Output a list of all successfully downloaded samples we have so far.
if false
    load NCBI_SRA_Mmusculus_download_04June2015_prelim.mat
    d = dataset('file', 'successfully_processed_list_11July2015_download.txt', 'ReadObsNames', false, 'ReadVarNames', false);
    
    ids1 = s.ids';
    ids2 = d.Var1;
    
    dlmcell('all_successfully_downloaded_up_through_11July2015_download.txt', [ids1;ids2]);
    
    return
end


% Pre-process, quality filter, and quality check the data.
if false
    qfparams.MRTHRESH = 0.7;
    qfparams.RCTHRESH = 4e6;
    qfparams.CORRCUTOFF = 0.55;
    qfparams.NZCUTOFF = 0.3; 
    qfparams.TPMCUTOFF = 1;
    
    if true
        load('NCBI_SRA_Mmusculus_full_data_up_to_19Sept2015.mat');
        [Y,sids,tids] = quality_filter(s, qfparams, 'Mmusculus');
    end
    
    if false
        plot_quality_filter_report('Mmusculus_quality_filter_summary_data.mat', qfparams, 'Mmusculus');
    end
    
    % Prepare the query table
    qt_full = read_ncbi_sra_query_table(qtfile);
    
    % Some entries from the SRA have been removed before we updated
    % the query table. Remove these samples from analysis.
    k = steq(sids, get(qt_full, 'ObsNames'));
    Y = Y(:,k);
    sids = sids(k);
    
    qt = qt_full(sids,:);
    save(mainDataFile, 'Y', 'tids', 'sids', 'qt', '-v7.3');
    
    return;
end



load(mainDataFile);
if false
    lY = lag_dataset(Y'.*repmat(qt.spots,1,size(Y,1))/1000000, qt.spots/1000000);
    save(mainDataFile, 'lY', '-append');
    return;
end

qt = NCBI_SRA_Mmusculus_build_and_analyze_query_table( qtfile ); % Load the query table
qt = qt(sids,:);

% reformat transcript IDs.
splits = regexpi(tids, ',', 'split');
sp = reshape([splits{:}], 2, length(tids))';
stids = tids;
tids = sp(:,2);

%lY = log10(Y' + 0.1);


% Number of HQ transcriptomes in the SRA as a function of time.
if false
    plot_ncbi_sra_growth(qt);
    plotSave('figures/growth_over_time/NCBI_SRA_Mmusculus_growth_v_time.png');
    close
end

% PCA. Percent variation explained vs. eigengene.
if false;
    [coef, pexp] = pexp_vs_components(lY, 'Mmusculus'); 
end

% PCA first 1-3 dimensions.
if false
    load('NCBI_SRA_Mmusculus_PCA_pexp_vs_eigengene_params.mat');
    NCBI_SRA_Mmusculus_plot_PCA( lY, coef, qt, pexp )
    s = lY*coef(:,1:100);
    pc_convergence_plot( s, qt, pexp )
end


% PCA plots by submission, for supplemental figures
if false;
    load('NCBI_SRA_Mmusculus_PCA_pexp_vs_eigengene_params.mat');
    pca_by_submission(lY, qt, coef, pexp);
    plotSave('figures/pca/NCBI_SRA_Mmusculus_PCA_by_submission.png');
    close;
end
    

% Marker OMP decomposition on FULL data.
% Run on $pw
if false
    punexp = 0;
    maxfeats = 500;
    somp = marker_OMP(standardize(lY), punexp, 'savememory', true, 'maxfeatures', maxfeats);
    save(sprintf('NCBI_SRA_Mmusculus_marker_OMP_decomposition_punexp_%0.2f_maxfeats_%0.0f.mat', punexp, maxfeats), 'somp');
end

% Compress the FULL data.
% Used for figure 2.
if false
    NMARKERS = 100;
    load('NCBI_SRA_Mmusculus_marker_OMP_decomposition_punexp_0.00_maxfeats_500.mat');
    model = tratrain(lY, lY(:, somp.S(1:NMARKERS)));
    model.reconstruction = lY(:, somp.S(1:NMARKERS))*model.b + repmat(model.b0,size(lY,1),1);
    save(['NCBI_SRA_Mmusculus_compression_and_reconstruction_nmarkers_', num2str(NMARKERS), '.mat'], 'model')
end


% Heatmap of original vs. reconstruction of FULL data.
% Part of Figure 2.
if false
    NMARKERS = 100;
    load('NCBI_SRA_Mmusculus_marker_OMP_decomposition_punexp_0.00_maxfeats_500.mat');
    load(['NCBI_SRA_Mmusculus_compression_and_reconstruction_nmarkers_', num2str(NMARKERS), '.mat'])
    
    model.reconstruction = tradict( lY(:,somp.S(1:NMARKERS)), model);

    heatmap_raw_vs_reconstructed(Y',somp, model, 'Mmusculus', false, NMARKERS);
    
    f = pred_v_actual_density_plot(lY, model.reconstruction, 'subsample_genes', 5000);
    hcb=colorbar;
    ca = caxis;
    set(hcb,'YTick',round(100*linspace(min(ca), max(ca), 4))/100);
    sf = get_standard_figure_font_sizes;
    set(gca, 'FontSize', sf.axis_tick_labels);
    plotSave('figures/heatmap_original_vs_reconstruction/insample_density_plot.png');
    
end


% Gene set analysis.
if true
    load('~/GitHub/transcriptome_compression/analysis/gene_ontology/Mmusculus_representative_gene_set_02-Apr-2016.mat');

    if true
        if false
            t = (Y').*repmat(qt.spots/1000000,1, size(Y,1) );
            o = qt.spots/1000000;
            nfolds = 20;

            results = evaluate_prospective_performance_3( t, o, tids, sets, qt, nfolds );
            save('NCBI_SRA_Mmusculus_evaluate_prospective_performance_3_results.mat', 'results', '-v7.3');
        else
            load('NCBI_SRA_Mmusculus_evaluate_prospective_performance_3_results.mat');
            analyze_prosperf3_results( results, qt );
        end
    end
    
    % Training on the full data.
    % This is the full model to use for future experiments.
    if false
        l = logspace(-6,-2,19);
        params = {l(7), false, false}; %l(7) selected from earlier cross validation.
        trainfun = @ridgefit_train;
        predfun = @ridgefit_predict;
        nmarkers = 100;
        
        % Test expression optimization.
        if false
            meanexp = mean(lY);
            [sY, model.train_mu, model.train_sig] = standardize(lY);
            model = geneset_cluster( sY, tids, sets, 'stats', model );
            model_expopt_100 = geneset_encode(sY, nmarkers, model, 'expression_delta', 100, 'mean_expression', meanexp);
            model_expopt_20 = geneset_encode(sY, nmarkers, model, 'expression_delta', 20, 'mean_expression', meanexp);
            model_expunopt = geneset_encode(sY, nmarkers, model, 'expression_delta', 0, 'mean_expression', meanexp);
            
            
            
            figure
            hold on
            plot(model_expopt_100.punexp, '-b')
            plot(model_expopt_20.punexp, '-r')
            plot(model_expunopt.punexp, '-k')

            
            figure;
            subplot(1,3,1);
            bar( 10.^mean(lY(:,model_expunopt.S)) - 0.1);
            axis square
            
            subplot(1,3,2);
            bar( 10.^mean(lY(:,model_expopt_20.S)) - 0.1);
            axis square
            
            subplot(1,3,3);
            bar(10.^mean(lY(:,model_expopt_100.S)) - 0.1);
            axis square
            
        end
        
        % expression optimization code suggests a delta (neighborhood size) of 20 is best.
        model = tradict_train( lY, tids, sets, trainfun, params, 'expression_delta', 20 );
        
        save('NCBI_SRA_Mmusculus_final_tradict_model.mat', 'model');
    end
    
    
    % Full cross validation.
    if false 
        trainfun = @ridgefit_train;
        predfun = @ridgefit_predict;
        l = logspace(-6,-2,19);
        params = {l(7), false, false};
        nfolds = 100;
        results = evaluate_prospective_performance_2( ... 
            lY, tids, sets, qt, predfun, trainfun, params, nfolds );
        
        save('NCBI_SRA_Mmusculus_evaluate_prospective_performance_2_results.mat', 'results');
    end
    
    
    % Results of full cross validation.
    if false
        NSTHRESH = 6;
        [ts, ns] = subadjust(results.target_proc, qt.Submission);
        ps = subadjust(results.pred_proc, qt.Submission);
        
        % gene set median relative errors. 
        rel = abs((results.pred_proc - results.target_proc)./results.target_proc);
        rel(rel == inf) = nan;
        mrel = nanmedian(rel); 
        
        % gene set Rsq between standardized expression.
        rsq = rsq_and_slope(ts(ns>=NSTHRESH,:),ps(ns>=NSTHRESH,:));
        
        % Gene set prediction error as a function of gene set size and
        % average expression.
        if true
            report = cell(length(sets),6);
            for i = 1 : length(sets)
                mask = steq(tids, sets{i});
                me = mean(lY(:,mask));
                avgexp = mean(me);
                seexp = std(me)/sqrt(sum(mask)-1);
                ngenes = sum(mask);

                report{i,1} = setnames{i,1};
                report{i,2} = ngenes;
                report{i,3} = avgexp; %sprintf('%0.2f +/- %0.2f', avgexp, seexp);
                report{i,4} = seexp;
                report{i,5} = mrel(i);
                report{i,6} = rsq(i);
            
                
            end
            
            dout = cell2dataset(report, 'VarNames', {'Set_Name', ...
                'Set_Size', 'Avg_Member_Exp', ...
                'StdErr_Member_Exp', 'Median_Relative_Error', ...
                'Intra_Submission_PCC'});
            
            
            export(dout, 'file', 'Mmusculus_geneset_accuracy_report.txt');
        end
        
    end
    
    
    
    if false % DEPRECATED
        rng('default')
        [sY, stats.train_mu, stats.train_sig] = standardize(lY);
        stats = geneset_cluster( sY, tids, sets, 'stats', stats );
        stats = geneset_encode(sY, 100, stats);

        save('NCBI_SRA_Mmusculus_geneset_encoding.mat', 'stats');
    end
    
    
    % Ridge fit training for predicting gene sets.
    if false % DEPRECATED
        load('NCBI_SRA_Mmusculus_geneset_encoding.mat');
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
        
    
        save('NCBI_SRA_Mmusculus_geneset_encoding.mat', 'stats', 'cvs');
    end
    
    
    if false % DEPRECATED
        load('NCBI_SRA_Mmusculus_geneset_encoding.mat');
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
        evaluate_prospective_performance(lY,qt, 'NCBI_SRA_Mmusculus_prospective_performance.mat');
    else
        load('NCBI_SRA_Mmusculus_prospective_performance.mat');
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

% Context specific performance. 
if false
    label = 'hematopoetic_lymphatic';
    class = {'hematopoetic', 'lymphatic'};
    %label = 'nervous';
    %class = {'nervous', 'developing_nervous'};
    
    load('NCBI_SRA_Mmusculus_PCA_pexp_vs_eigengene_params.mat');
    saveFile = ['NCBI_SRA_Mmusculus_context_specific_prospective_performance_', label, '.mat'];
    NCBI_SRA_Mmusculus_context_specific_performance_trial( lY, qt, coef, class, saveFile );
    plotSave(['figures/prospective_performance/NCBI_SRA_Mmusculus_context_specific_prospective_performance_', label, '.png']);
    close
    
%     if false
%         load('Mmusculus_context_specific_test_data.mat');
%         pred_v_actual_density_plot(ytrues,yhats);
%         plotSave('figures/prospective_performance/test_context_specific_hematopoetic.png');
%         iminvert('figures/prospective_performance/test_context_specific_hematopoetic.png');
%         close
%         
%         ytrues = [randn(100, 400); 2.5+randn(200, 400)];
%         yhats = ytrues + 1.2*randn(300,400);
%         pred_v_actual_density_plot(ytrues,yhats);
%        
%         plotSave('figures/prospective_performance/test_context_specific_nervous.png');
%         iminvert('figures/prospective_performance/test_context_specific_nervous.png');
%         close
%         
%     end
    
end

t = toc;
if completiontext
    fprintf('Sending text.\n');
    msg = sprintf('Job took %0.2f seconds.\n... you''re awesome.', t);
    send_text_message('919-757-1609', 'AT&T', 'Job Finished', msg)
end


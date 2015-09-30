clear;
rng('default');

% Initialize environment
repo = 'transcriptome_compression/';
datadir = ['~/GitHub/data/', repo, 'Mmusculus/'];
path(genpath(['~/GitHub/', repo]), path);
cd(datadir);

% General variables
qtfile = 'Mmusculus_query_table_04June2015_.csv';
mainDataFile = 'NCBI_SRA_Mmusculus_download_04June2015_prelim_processed.mat';

% Construction full collection of raw transcriptomes.
if true
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
    load NCBI_SRA_Mmusculus_download_04June2015_prelim.mat
    NCBI_SRA_Mmusculus_preprocess(s,qt,mainDataFile);
    return;
end

load(mainDataFile);
qt = NCBI_SRA_Mmusculus_build_and_analyze_query_table( qtfile ); % Load the query table
qt = qt(sids,:);

lY = log10(Y' + 0.1);


% Number of HQ transcriptomes in the SRA as a function of time.
if false
    plot_ncbi_sra_growth(qt);
    plotSave('figures/growth_over_time/NCBI_SRA_Mmusculus_growth_v_time.png');
    close
end

% PCA. Percent variation explained vs. eigengene.
if false;
    load('NCBI_SRA_Mmusculus_PCA_pexp_vs_eigengene_params.mat');
    [coef, pexp] = pexp_vs_components(lY, 'Mmusculus'); 
end

% PCA first 1-3 dimensions.
if false
    load('NCBI_SRA_Mmusculus_PCA_pexp_vs_eigengene_params.mat');
    NCBI_SRA_Mmusculus_plot_PCA( lY, coef, qt, pexp )
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

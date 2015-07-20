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

% Load the query table
if true; % Leave true
    qt = NCBI_SRA_Mmusculus_build_and_analyze_query_table( qtfile );
end

% Pre-process, quality filter, and quality check the data.
if false
    load NCBI_SRA_Mmusculus_download_04June2015_prelim.mat
    NCBI_SRA_Mmusculus_preprocess(s,qt,mainDataFile);
    return;
end

load(mainDataFile);
lY = log10(Y' + 0.1);

% PCA. Percent variation explained vs. eigengene.
if false;
    NCBI_SRA_Mmusculus_pexp_vs_components(lY);
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
    heatmap_raw_vs_reconstructed(Y',somp, model, 'Mmusculus', false, NMARKERS);
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


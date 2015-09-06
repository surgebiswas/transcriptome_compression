clear;
rng('default');

% Initialize environment
repo = 'transcriptome_compression/';
datadir = ['~/GitHub/data/', repo, 'Athaliana/'];
path(genpath(['~/GitHub/', repo]), path);
cd(datadir);

mainDataFile = 'NCBI_SRA_Athaliana_full_data_up_to_18May2015_processed_updated_09June2015.mat';
queryTable = 'Athaliana_query_table_18May2015_shuffled_.csv';

load(mainDataFile);
lY = log10(Y' + 0.1);

load('NCBI_SRA_Athaliana_PCA_pexp_vs_eigengene_params.mat');
load(sprintf('NCBI_SRA_%s_marker_OMP_decomposition_punexp_0.00_maxfeats_500.mat', 'Athaliana'));

% First let's see how the markers we select using the PC space compares to
% the markers we select using the full expression matrix.
s = lY*coef(:,1:500);

maxfeats = 100;
pstats = pc_omp(lY, s, 'maxfeats', maxfeats);

lambda = logspace(-8, 1, 20);
mpc = tratrain(lY, lY(:,pstats.L), 'lambda', lambda);
mfull = tratrain(lY, lY(:,somp.S(1:maxfeats)), 'lambda', lambda);
mrand = tratrain(lY, lY(:,randsample(size(lY,2),maxfeats)), 'lambda', lambda);

notBoxPlot([mrand.mse(:,mrand.lambda_star_idx), mpc.mse(:,mpc.lambda_star_idx), mfull.mse(:,mfull.lambda_star_idx)])


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
if true; 
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
if true;
    NCBI_SRA_Mmusculus_pexp_vs_components(lY);
end
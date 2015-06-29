clear;
rng('default');

% Initialize environment
repo = 'transcriptome_compression/';
datadir = ['~/GitHub/data/', repo, 'Mmusculus/'];
path(genpath(['~/GitHub/', repo]), path);
cd(datadir);

% General variables
qtfile = 'Mmusculus_query_table_04June2015_.csv';

% Prepare and analyze the query table
if true; 
    qt = NCBI_SRA_Mmusculus_build_and_analyze_query_table( qtfile );
end


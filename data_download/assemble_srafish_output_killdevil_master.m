path(genpath('/proj/dangl_lab/sbiswas/GitHub/transcriptome_compression/'), path);
saveFile = '/proj/dangl_lab/sbiswas/GitHub/data/transcriptome_compression_server_only/Mmusculus/download_11July2015/Mmusculus_download_11July2015_assembled_srafish_output.mat';

% Successfully processed list of SRA runs.
% Should have a variable named 'sp' for successfully processed IDs. This
% should be a cell aray of SRA IDs. 
load('/proj/dangl_lab/sbiswas/GitHub/data/transcriptome_compression_server_only/Mmusculus/download_11July2015/successfully_processed_list_11July2015_download.mat');

% Location of the output directory where all SRA ID output directories are
% kept.
srafish_outdir = '/lustre/scr/s/b/sbiswas/NCBI_SRA/Mmusculus/';

s = assemble_srafish_output(srafish_outdir, sp);

save(saveFile, 's', '-v7.3');



clear;
old = cd('~/Documents/surge/science/gene_ontology');

GO = geneont('live',true, 'tofile', ['gene_ontology_db_', datestr(today), '.obo']); % this step takes a while
cd(old);

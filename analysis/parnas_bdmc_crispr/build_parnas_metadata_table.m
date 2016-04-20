% Builds metadata table for Parnas et al 2015 RNA seq data.

clear;
rng('default');
cd('/home/sbiswas/GitHub/data/transcriptome_compression/parnas_bdmc_crispr');

q = read_ncbi_sra_query_table('SRA248232_query_table.csv');

rawmd = cell(size(q,1),5);
for i = 1 : size(q,1)
    g = get_geo_metadata_parnas(q.SampleName{i});
    rawmd{i,1} = g.ko;
    rawmd{i,2} = g.tcl;
    rawmd{i,3} = g.con;
    rawmd{i,4} = g.treat;
    rawmd{i,5} = g.source_name;
    disp(i)
end

% parse out source name some more to collect what looks like
% the library name and potentially plate ID??
lib = cell(size(q,1),1);
unknown1 = cell(size(q,1),1);
for i = 1 : size(q,1)
    disp(i)
    m = regexpi(rawmd{i,5}, '(MGLib.)', 'tokens');
    if isempty(m)
        lib{i} = 'NA';
    else
        lib{i} = m{1}{1};
    end
    
    
    s = rawmd{i,5};
    s = regexprep(s, rawmd{i,1}, '');
    s = regexprep(s, '(MGLib.)', '');
    s = regexprep(s, '\.\.', '');
    s = regexprep(s, '\_', '');
    s = regexprep(s, 'T\d$', '');

    unknown1{i} = s;
end


xdu = cell2dataset([rawmd, lib, unknown1], 'VarNames', {'KO', 'line', 'condition', 'treatment', 'source_name', 'lib', 'unknown'}, ...
    'ObsNames', get(q, 'ObsNames'));

save('parnas_design.mat', 'xdu');
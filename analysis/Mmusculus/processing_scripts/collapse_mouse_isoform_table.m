function [ctids, cM, gids, comgids] = collapse_mouse_isoform_table(tids, M)
% tids = transcript ID's
% M     = transcripts x samples expression matrix
% 
% Assumes the gencode M5 annotation.
% Example transcript ID:
% ENSMUST00000194992.2|ENSMUSG00000025900.7|OTTMUSG00000049985.2|OTTMUST00000127194.1|Rp1-002|Rp1|3047|UTR5:1-54|CDS:55-912|UTR3:913-3047|
%
% citds = collapsed transcript ID's. Retains only the Gene ID and the gene
% common name.
% cM = collapsed expression table.
% gids = extracted gene IDs for each transcript (length == length(tids));
% comgids = extracted gene common names for each transcript ''.


gids = cell(length(tids),1);
comgids = cell(length(tids),1);
splits = regexpi(tids, '\|', 'split');
for i = 1 : length(gids)
    gids{i} = splits{i}{2};
    comgids{i} = splits{i}{6};
end

[ug, ~, ic] = unique(gids);
cM = zeros(length(ug),size(M,2));
ctids = cell(length(ug),1);
for i = 1 : length(ug)
    cM(i,:) = sum(M(ic == i, :),1);
    
    ctids{i} = strjoin([ug(i); unique(comgids(ic == i))]', ',');    
    
    if mod(i,1000) == 0
        fprintf('%0.0f%% collapsed.\n', 100*(i/length(ug)));
    end
end
fprintf('%0.0f%% collapsed.\n', 100*(i/length(ug)));




end
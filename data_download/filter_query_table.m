function qt = filter_query_table( query_table, varargin )
% query_table = query table to filter.
%
% processed_list = text file of a list of samples that have been processed.
% e.g 
% DRR008476
% DRR008477
% DRR008478
% ERR229826
% ERR229827
% ...
%

minreads = setParam(varargin, 'minreads', 0);
platform = setParam(varargin, 'platform', []);
processed_list = setParam(varargin, 'processed_list', []);

if ~isempty(processed_list) 
    d = dataset('file', processed_list, 'ReadObsNames', false, 'ReadVarNames', false);
    processed_list = d.Var1;
    
    plhash = containers.Map;
    for i = 1 : length(processed_list)
        plhash(processed_list{i}) = 1;
    end
end




qt = read_ncbi_sra_query_table(query_table);
oqt = get(qt, 'ObsNames');
tokeep = true(size(qt,1),1);
for i = 1 : size(qt,1)
    if qt.spots(i) < minreads
        tokeep(i) = false;
    end
    
    if ~isempty(platform) && ~strcmpi(qt.Platform{i}, platform)
        tokeep(i) = false;
    end
    
    if ~isempty(processed_list) && plhash.isKey(oqt{i})
        tokeep(i) = false;
    end
end

fid = fopen(query_table);
fout = fopen([query_table, '.filtered'], 'w+');

l = fgetl(fid); % Header line
fprintf(fout, '%s\n', l); % Header line

lines = cell(size(qt,1),1);
for i = 1 : length(lines)
    lines{i} = fgetl(fid);
end
fclose(fid);

kept_lines = lines(tokeep);
for i = 1 : length(kept_lines)
    fprintf(fout, '%s\n', kept_lines{i});
end




end


function [pathway_data, sets, setnames] = build_kegg_pathway_gene_sets( organism )
% http://www.mathworks.com/help/bioinfo/examples/connecting-to-the-kegg-api-web-service.html

switch organism
    case 'Mus musculus'
        orgcode = 'mmu';
    case 'Arabidopsis thaliana'
        orgcode = 'ath';
end

base = 'http://rest.kegg.jp/';

% Get list of pathways. These are the set names.
operation = 'list/';
database = 'pathway/';
organismCode = orgcode;
pathway_list = urlread(strcat(base,operation,database,organismCode));
pathway_list = regexpi(pathway_list,'[^\n]+','match')'; % convert to cellstr
num_pathways = numel(pathway_list); % total number of pathways

fprintf('%0.0f pathways found.\n', num_pathways);


% Now get the list of genes associated with each pathway.
sets = cell(length(pathway_list),1);
setnames = sets;
for i = 1 : length(pathway_list)
    fprintf('%0.0f. Retrieving data for pathway: %s\n', i, pathway_list{i});
    
    operation = 'get/';
    pid = regexpi(pathway_list(i),'(?<=path:)\w+','match');
    record = urlread(char(strcat(base,operation,pid{1})));
    
    % Retrieve the KO_PATHWAY id and all other alias pathway entries for it.
    ko_id = regexpi(record,'(?<=KO\w+PATHWAY\s+)\w*','match');
    operation = 'link/';
    database = 'pathway/';
    allPathwayIDs = urlread(strcat(base,operation,database,ko_id{1}));
    map_id = regexpi(allPathwayIDs,'(?<=\w+\W+\w+\s+path:)(?=map)\w*', 'match');
    
    % Get gene list
    try
        [gene_common_ids, gene_common_names] = parse_gene_names(record);
        sets{i} = gene_common_ids;
        setnames{i} = gene_common_names;

        pathway_data(i).name = pathway_list{i};
        pathway_data(i).gene_common_ids = gene_common_ids;
        pathway_data(i).gene_common_names = gene_common_names;
        pathway_data(i).map_id = map_id;
    catch
        warning(sprintf('Could not parse data for pathway'));
        disp(record);
    end
end
torm = cellfun(@isempty, sets);
sets(torm) = [];
setnames(torm) = [];
pathway_data(torm) = [];



    function [gene_common_ids, gene_common_names] = parse_gene_names(record)
        f1 = regexprep(record, '.+?GENE', '');
        
        q = char(regexp(f1, '\n', 'split')');
        ci = find(q(:,1) ~= ' ', 1);
        
        f1(ci:end,:) = [];
                

        splits = strtrim(regexp(f1, '\n', 'split'))';
        splits(cellfun(@isempty, splits)) = [];

        gc = regexp(splits, '\d\s+(.+?);', 'tokens');
        ggc = [gc{:}]';
        gene_common_ids = [ggc{:}]';


        gc = regexp(splits, ';(.+)', 'tokens');
        ggc = [gc{:}]';
        gene_common_names = strtrim([ggc{:}]');
    end


end


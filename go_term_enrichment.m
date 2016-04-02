function report = go_term_enrichment(genes, short_list, organism, varargin)
% aspect = 'C', 'F', or 'P'.
% organism = 'Arabidopsis thaliana', or 'Mus musculus'
% genes = all gene IDs
% short_list = list of gene IDs of interest (e.g. those in a cluster).

    %%% VARIABLE ARGUMENTS %%% 
    GO = setParam(varargin, 'GOdb', []);
    pvalue_cutoff = setParam(varargin, 'pvalue_cutoff', 0.01); % Terms to include in report.
    
    % Unlike 'amap', this hash maps a gene ID to not only its go terms but
    % also all of its relatives. This can be expensive to compute for many
    % go terms, so a precomputed version can be supplied 
    relhash = setParam(varargin, 'relatives_hash', []);
    
    % Set up
    if isempty(GO)
        fprintf('Loading GO database ... ');
        GOdbpath = '/Users/sbiswas/Documents/surge/science/gene_ontology/gene_ontology_db_17-Feb-2016.obo';
        GO = geneont('File', GOdbpath);
        fprintf('Done.\n');
    end
    
    switch organism
        case 'Arabidopsis thaliana'
            matfile = '/Users/sbiswas/Documents/surge/science/gene_ontology/Athaliana_gene_association_16Feb2016.mat';
        case 'Mus musculus'
            matfile = '/Users/sbiswas/Documents/surge/science/gene_ontology/Mmusculus_gene_association_29Mar2016.mat';
    end
    load(matfile); % see build_association_mat_file.m
    
    
    % Perform enrichment analysis
    
    % Build relatives hash or use supplied.
    if isempty(relhash)
        fprintf('Computing relatives hash ... ');
        relhash = containers.Map;
        for i = 1 : numel(genes)
            if isKey(amap, genes{i})
                relhash(genes{i}) = getrelatives(GO,amap(genes{i})); % slow step, rate limiting.
            end
            if mod(i, round(numel(genes)/10)) == 0
                fprintf('%0.2f%% ... ', 100*i/numel(genes));
            end
        end
        fprintf('Done.\n');
    end
    
    
    
    % Perform counting.
    fprintf('Performing counting ... ');
    m = GO.Terms(end).id;           % gets the last term id
    totalcount = zeros(m,1);        % a vector of GO term counts for all genes.
    shortlistcount = zeros(m,1);    % a vector of GO term counts for interesting genes.
    for i = 1:numel(genes)
        if isKey(relhash,genes{i})
            goid = relhash(genes{i}); %getrelatives(GO,amap(genes{i})); 
            
            % update vector counts
            totalcount(goid) = totalcount(goid) + 1;
            if any( strcmpi(genes{i}, short_list) )
               shortlistcount(goid) = shortlistcount(goid) + 1;
            end
        end
        
        
    end
    fprintf('Done.\n');
    
    % Perform significance testing
    fprintf('Performing significance testing ... ');
    pvalues = go_term_enrichment_pvalues(shortlistcount, totalcount);
    fprintf('Done.\n');
    
    % Generate report
    fprintf('Generating report ... ');
    [pvs,idx] = sort(pvalues);
    pvs(pvs > pvalue_cutoff) = [];
    idx(pvs > pvalue_cutoff) = [];
    
    rep = cell(length(pvs), 4);
    for i = 1 : length(pvs)
        term = idx(i);
        rep{i,1} = char(num2goid(term));
        rep{i,2} = pvs(i);
        rep{i,3} = sprintf('(%0.0f/%0.0f)/(%0.0f/%0.0f)', shortlistcount(term), numel(short_list), totalcount(term), numel(genes));
        rep{i,4} = GO(term).term.name;
    end
    
    report = cell2dataset(rep, 'VarNames', {'GO_Term', 'p_value', 'enrichment', 'description'});
    
    fprintf('Done.\n');
    
end


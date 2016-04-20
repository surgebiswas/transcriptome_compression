function model = tradict_train( lY, tids, sets, trainfun, trainparams, varargin )

    nmarkers = setParam(varargin, 'nmarkers', 100);
    expdelta = setParam(varargin, 'expression_delta', 0);
    
    
    
    % perform the endcoding
    meanexp = mean(lY);
    [sY, model.train_mu, model.train_sig] = standardize(lY);
    model = geneset_cluster( sY, tids, sets, 'stats', model );
    model = geneset_encode(sY, nmarkers, model, 'expression_delta', expdelta, 'mean_expression', meanexp);
    
    
    markers = model.S;
    target = model.geneset.sy_sets;
    
    model.fit = trainfun(lY(:,markers), target, trainparams);
    model.fit_allgenes = trainfun(lY(:,markers), lY, trainparams);
    
end


function model = tradict_train( lY, tids, sets, trainfun, trainparams, varargin )

    nmarkers = setParam(varargin, 'nmarkers', 100);
    
    
    % perform the endcoding
    [sY, model.train_mu, model.train_sig] = standardize(lY);
    model = geneset_cluster( sY, tids, sets, 'stats', model );
    model = geneset_encode(sY, nmarkers, model);
    
    
    markers = model.S;
    target = model.geneset.sy_sets;
    
    model.fit = trainfun(lY(:,markers), target, trainparams);
    
end


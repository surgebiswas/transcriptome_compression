function model = tradict_train_pmvn( t, o, tids, sets, varargin )

    nmarkers = setParam(varargin, 'nmarkers', 100);
    expdelta = setParam(varargin, 'expression_delta', 0);

    % Perform lag
    [zlag, model.lag_priors] = lag_dataset(t, o);
    
    % perform the endcoding
    meanexp = mean(zlag);
    [sY, model.train_mu, model.train_sig] = standardize(zlag);
    model = geneset_cluster( sY, tids, sets, 'stats', model );
    model = geneset_encode(sY, nmarkers, model, 'expression_delta', expdelta, 'mean_expression', meanexp);
    markers = model.S;
    
    % Learn z_m, \mu^{(m)}, and \Sigma^{(m)}
    [model.fit.markers.z, model.fit.markers.mu, model.fit.markers.Sigma] = ...
            learn_pmvn(t(:,markers), o, zlag(:,markers));
    zlag(:,markers) = model.fit.markers.z; % update lag estimates for markers.

    % Learn mean and cross covariances for gene sets (processes)
    model.fit.geneset.mu = mean(model.geneset.sy_sets);
    model.fit.geneset.cross_cov = cross_cov(model.fit.markers.z, model.geneset.sy_sets);
    
    % Learn mean and cross covariances for all genes
    model.fit.genes.mu = mean(zlag);
    model.fit.genes.var = var(zlag);
    model.fit.genes.cross_cov = cross_cov(model.fit.markers.z, zlag);
    
    function c = cross_cov(x,y)
        c = bsxfun(@minus,x,mean(x))'*bsxfun(@minus,y,mean(y))/(size(x,1)-1);
    end

end


function compression_ensemble_train( lY, varargin )

% if a save file is specified, progress is saved for quick restarts.
saveFile = setParam(varargin, 'savefile', []); 
numpcacomponents = setParam(varargin, 'numpcacomponents', 10);
minclusterdepth = setParam(varargin, 'minclusterdepth', 100);
nclusterboot = setParam(varargin, 'nclusterboot', 100);
subresidual = setParam(varargin, 'subresidual', 0.01);
numanchors = setParam(varargin, 'numanchors', 5);
nummarkers = setParam(varargin, 'nummarkers', 500);


% Check progress.
compute_pca = true;
compute_clusters = true;
compute_anchors = true;
compute_csmarkers = true;
if ~isempty(saveFile)
    if exist(saveFile, 'file')
        load(saveFile)
        compute_pca = ~exist('s', 'var');
        compute_clusters = ~exist('cidx', 'var');
        compute_anchors = ~exist('anchor_stats', 'var');
        compute_csmarkers = ~exist('cs_stats', 'var');
    end
end

% 1. PCA computation to define a low dimensional expression space.
if compute_pca
    fprintf('Computing PCA.\n');
    [coef, ~, ~, ~, pexp] = pca(lY, 'NumComponents', numpcacomponents);
    s = lY*coef; % PC scores
    save(saveFile, 's', 'coef', 'pexp');
end

% 2. Partition the expression space into clusters.
if compute_clusters
    fprintf('Computing expression clusters.\n');
    cidx = pinch_clusters(s, kmeans_conclust(s, 'nrep', nclusterboot), minclusterdepth);
    save(saveFile, 'cidx', '-append');
end

% 3. Define anchor markers
if compute_anchors
    fprintf('Computing anchor markers.\n');
    sY = standardize(lY);
    anchor_stats = marker_OMP(sY, 0, 'maxfeatures', numanchors, 'subresidual', subresidual);
    anchor_beta = pinv(sY(:,anchor_stats.S))*sY;
    
    save(saveFile, 'anchor_stats', 'anchor_beta', '-append');
end

% 4. Define context-specific markers.
if compute_csmarkers
    fprintf('Computing context-specific markers.\n');
    
    R = sY - sY(:,anchor_stats.S)*anchor_beta; % Remaining residual for context specific decomposition.
    cs_stats = cell(max(cidx),1);
    budget = nummarkers - numanchors;
    for i = 1 : max(cidx)
        nmarkers = round(budget*sum(cidx==i)/length(cidx));
        
        % Perform context-specific decompositions on the residual after
        % computing the anchors.
        cs_stats{i} = marker_OMP(R(cidx==i,:), 0, 'maxfeatures', nmarkers, 'subresidual', 0.01);
    end
    save(saveFile, 'cs_stats', '-append');
end

% 5. Learn the kernel covariance 



end


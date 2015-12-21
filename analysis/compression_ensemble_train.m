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
compute_kernel = true;
if ~isempty(saveFile)
    if exist(saveFile, 'file')
        load(saveFile)
        compute_pca = ~exist('s', 'var');
        compute_clusters = ~exist('cidx', 'var');
        compute_anchors = ~exist('anchor_stats', 'var');
        compute_csmarkers = ~exist('cs_stats', 'var');
        compute_kernel = ~exist('model', 'var');
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
[sY, mu, sig] = standardize(lY);
if compute_anchors
    fprintf('Computing anchor markers.\n');
    anchor_stats = marker_OMP(sY, 0, 'maxfeatures', numanchors, 'subresidual', subresidual);
    anchor_beta = pinv(sY(:,anchor_stats.S))*sY;
    
    save(saveFile, 'anchor_stats', 'anchor_beta', 'mu', 'sig', '-append');
end

% 4. Define context-specific markers.
R = sY - sY(:,anchor_stats.S)*anchor_beta; % Remaining residual for context specific decomposition.
if compute_csmarkers
    fprintf('Computing context-specific markers.\n');
    
    cs_stats = cell(max(cidx),1);
    budget = nummarkers - numanchors;
    for i = 1 : max(cidx)
        nmarkers = round(budget*sum(cidx==i)/length(cidx));
        
        % Perform context-specific decompositions on the residual after
        % computing the anchors.
        cs_stats{i} = marker_OMP(R(cidx==i,:), 0, 'maxfeatures', nmarkers, 'subresidual', 0.01, 'storecoefficients', true);
    end
    save(saveFile, 'cs_stats', '-append');
end

% 5. Learn the kernel covariance 
if compute_kernel
    model.mu = mu;
    model.sig = sig;
    model.b_anchor = anchor_beta;
    model.idx_anchor = anchor_stats.S;

    model.idx_all = model.idx_anchor;
    for i = 1 : length(cs_stats)
        model.idx_all = [model.idx_all, cs_stats{i}.S];
    end

    model.means = zeros(max(cidx), length(model.idx_all));
    for i = 1 : length(cs_stats)
        model.means(i,:) = mean(R(cidx==i,model.idx_all));
        model.b_cs{i} = cs_stats{i}.b;
        model.idx_cs{i} = cs_stats{i}.S;
    end

    %[model.alpha, model.alpha_scale] = learn_alpha(sY,model);
    model.alpha = initalpha(sY(:,model.idx_all));
    model.alpha_scale = 0.62; % learned from learn_alpha -- need to stream line this.
    
    
    save(saveFile, 'model', '-append'); 
end


    function [a,s] = learn_alpha(x, model)
        NFOLD = 10;
        NITR = 400;
        
        cvind = crossvalind('Kfold', size(x,1), NFOLD);
        alphas = cell(NFOLD,1);
        for k = 1 : max(cvind)
            alphas{k} = initalpha(x(cvind~=k,model.idx_all));
        end
        
        scales = zeros(1,NITR);
        s = 0.6198; %%%%%%
        testLogLike = -Inf;
        prevLL = testLogLike;
        a2 = 1;
        for itr = 1 : NITR
            if s > 0
                testLogLike = 0;
                for k = 1 : max(cvind)
                    % Training.
                    model.alpha = alphas{k};
                    model.alpha_scale = s;

                    % Test evaluation.
                    yhat = compression_ensemble_predict(x(cvind==k,model.idx_all), model);
                    ytrue = x(cvind==k,:);

                    testLogLike = testLogLike + -sum(sum( (yhat-ytrue).^2 ));
                end
            else
                testLogLike = -Inf;            
            end
            
            % Hastings update
            a1 = exp(testLogLike - prevLL);
            a = a1*a2; % Hastings ratio
            if a >= 1
                scales(itr) = s;
                prevLL = testLogLike;
            else
                if rand < a 
                    scales(itr) = s;
                else
                    scales(itr) = scales(itr-1);
                end
                prevLL = prevLL;
            end
            
            
            
            % New proposal
            s = 0.01*randn + scales(itr);
%             s = exprnd(scales(itr));
%             a2 = exppdf(scales(itr), s)/exppdf(s,scales(itr)); % Proposal probability ratio.

            fprintf('Itr: %0.0f\tScale: %0.4f\n', itr, scales(itr));
        end
        
        
    end

    function a = initalpha(x)
        a = inv(cov(x) + 0.05*eye(size(x,2))); % Need smarter regularization. 5% for now.
    end


end


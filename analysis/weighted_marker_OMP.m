function stats = weighted_marker_OMP( y, c, varargin )
% y = [samples x genes]

    savememory = setParam(varargin, 'savememory', false);
    maxfeatures = setParam(varargin, 'maxfeatures', Inf);
    storecrosscorr = setParam(varargin, 'storecrosscorr', false); % Store cross correlation between residual and original matrix?
    storecoefficients = setParam(varargin, 'storecoefficients', false); % Store OMP coefficients?
    w = setParam(varargin, 'obsweights', ones(size(y,1),1)/size(y,1));
    if ~iscolumn(w); w = w'; end
    
    if savememory; fprintf('Save memory? yes\n'); else fprintf('Save memory? no\n'); end
    fprintf('Max. features: %0.0f\n', maxfeatures);
    fprintf('Tolerated proportion of unexplained variance: %0.4f\n', c);
    
    
    
    sy = wst(y,w);
    wsy = sqrt(repmat(w,1,size(sy,2))).*sy;

    r = wsy;
    Phi = [];
    S = [];
    a = zeros(size(y));
    tvy = tv(wsy, w);
    punexp = tv(r, w)/tvy;
    crosscorr = [];
    while punexp(end) > c && length(S) < maxfeatures;
        
        % Naively computing the inner product will produce a genes x genes
        % correlation matrix, which may be quite difficult to store in
        % memory.
        if savememory
            [k,z] = max_sum_abs_inner_prod(r',wsy);
        else
            z = mean(abs(r'*wsy));
            [~, k] = max(z); 
        end
        
        if storecrosscorr; crosscorr = [crosscorr; z]; end;
        S = [S,k];
        Phi = [Phi, wsy(:,k)];
        
        b = Phi'*Phi \ Phi'*wsy;
        a = Phi*b;
        r = wsy - a;
        
        punexp(end+1) = tv(r,w)/tvy;
        
        fprintf('Iteration: %0.0f.\tProp. unexplained variance: %0.4f\n', length(punexp)-1, punexp(end));
    end

    stats.S = S;
    if storecoefficients; stats.b = b; else stats.b = []; end
    punexp(1) = []; % always == 1
    stats.punexp = punexp;
    stats.crosscorr = crosscorr;

    function [k, z] = max_sum_abs_inner_prod(rt,y)
        % The trick is to realize that we will sum over the rows of r'*y.
        % Thus the final vector is of dimension 1 x genes.
        % Thus we can compute the inner product in parts.
        BATCHSIZE = 100;
        z = zeros(1, size(y,2));
        
        ind = sort(crossvalind('Kfold', size(y,2), round(size(y,2)/100)));
        for j = 1 : max(ind)
            z(ind == j) = mean(abs(rt*y(:,ind ==j))); % mean taken across genes.
        end
        
        % Normalize inner products to be on the scale of correlations
        % z = z/size(y,2); % Dont need to normalize by size(y,1) since the matrices are weighted.
        [~,k] = max(z);
    end
    
    function v = tv(x,w)
        v = mean(var(x,w));
    end

    function sy = wst(y,w)
        % Weighted standardization.
        % sy has weighted mean = 0 and weighted stddev/variance = 1;
        sy = y;
        for i = 1 : size(y,2);
            sy(:,i) = (y(:,i) - w'*y(:,i))/sqrt(var(y(:,i),w));
        end
    end

end


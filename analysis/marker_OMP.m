function stats = marker_OMP( y, c, varargin )
% y = [samples x genes]

    savememory = setParam(varargin, 'savememory', false);
    maxfeatures = setParam(varargin, 'maxfeatures', Inf);
    storecrosscorr = setParam(varargin, 'storecrosscorr', false); % Store cross correlation between residual and original matrix?
    storecoefficients = setParam(varargin, 'storecoefficients', false); % Store OMP coefficients?
    subsampleresidual = setParam(varargin, 'subresidual', 0.1); % Fraction (0,1) of residual to subsample to. 
    pickrandomly = setParam(varargin, 'pickrandomly', false);
    
    if savememory; fprintf('Save memory? yes\n'); else fprintf('Save memory? no\n'); end
    fprintf('Max. features: %0.0f\n', maxfeatures);
    fprintf('Tolerated proportion of unexplained variance: %0.4f\n', c);
    
    
  
    r = y;
    Phi = [];
    S = [];
    a = zeros(size(y));
    tvy = tv(y);
    punexp = tv(r)/tvy;
    crosscorr = [];
    while punexp(end) > c && length(S) < maxfeatures;
        if pickrandomly
            k = randsample(setdiff(1:size(y,2), S), 1);
        else
        
            % Use r_eff? an effective, but reduced representation of the
            % residual.
            if true
                vr = var(r);
                %pct = prctile(vr, 90);
                ki = randsample(length(vr), round(length(vr)*subsampleresidual) );% vr > (pct - 1e-8); % Subtract 1e-8 to account for numerical precision errors in first iteration.

                r_eff = r(:,ki);
            else
                r_eff = r;
            end

            % Naively computing the inner product will produce a genes x genes
            % correlation matrix, which may be quite difficult to store in
            % memory.
            if savememory
                [k,z] = max_sum_abs_inner_prod(r_eff',y);
            else
                z = sum(abs(r_eff'*y))/(size(y,2)*size(y,1));
                [~, k] = max(z); 
            end

            if storecrosscorr; crosscorr = [crosscorr; z]; end;
        
        end
        S = [S, k];
        Phi = [Phi, y(:,k)];
        
        b = Phi'*Phi \ Phi'*y;
        a = Phi*b;
        r = y - a;
        
        punexp(end+1) = tv(r)/tvy;
        
        fprintf('Prop. unexplained variance: %0.4f\n', punexp(end));
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
            z(ind == j) = sum(abs(rt*y(:,ind ==j)));
        end
        % Normalize inner products to be on the scale of correlations
        z = z/(size(y,2)*size(y,1));
        [~,k] = max(z);
    end
    
    function v = tv(x)
        v = sum(var(x));
    end

end


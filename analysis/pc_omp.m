function stats = pc_omp( x, s, varargin )
% x = [samples x genes] expression matrix.
% s = [samples x principal components] scores matrix.

maxfeats = setParam(varargin, 'maxfeats', 100);
maxpropunexplained = setParam(varargin, 'punexp', 0);
w = setParam(varargin, 'obsweights', ones(size(x,1),1)/size(x,1));

assert(abs(sum(w) - 1) < 1e-10, 'Observation weights must sum to 1');

% Center and scale the gene expression matrix.
% Center the PC matrix. We only center the PC matrix because we want the
% first few components to retain their 'importance'.
x = wst(x,w,false);
s = wst(s,w,true);

L = [];
Phi = [];
tvs = tv(s,w);
r = s;
punexp = 1;
while punexp(end) > maxpropunexplained && length(L) < maxfeats;
    A = r'*weight(x,w); % A_{ij} = \sum_k r_{ki}*x_{ki}*w_k.
    z = mean(abs( A )); % Average absolute correlation.
    
    % Greedy selection.
    [~,l] = max(z);
    L = [L,l];
    Phi = [Phi, x(:,l)];
    
    % Orthogonal projection
    b = Phi'*Phi \ Phi'*s;
    a = Phi*b;
    r = s - a;
    
    punexp(end+1) = tv(r,w)/tvs;    
    fprintf('Iteration: %0.0f.\tProp. unexplained variance: %0.4f\n', length(punexp)-1, punexp(end));
    
end

stats.L = L;
stats.b = b;
punexp(1) = []; % always == 1
stats.punexp = punexp;



    function xw = weight(x,w)
        % weight rows of x.
        xw = repmat(w, 1, size(x,2)).*x;
    end

    function v = tv(x,w)
        v = mean(var(x,w),2);
    end

    function sy = wst(y,w,centeronly)
        % Weighted standardization.
        % sy has weighted mean = 0 and weighted stddev/variance = 1;
        sy = y;
        for i = 1 : size(y,2);
            if centeronly
                sy(:,i) = y(:,i) - w'*y(:,i);
            else
                sy(:,i) = (y(:,i) - w'*y(:,i))/sqrt(var(y(:,i),w));
            end
        end
    end



end


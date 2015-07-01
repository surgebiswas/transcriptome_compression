function stats = marker_OMP( y, c, varargin )
% y = [samples x genes]

    MAKEMOVIE = setParam(varargin, 'makemovie', false);
    savememory = setParam(varargin, 'savememory', false);
    maxfeatures = setParam(varargin, 'maxfeatures', Inf);
    storecrosscorr = setParam(varargin, 'storecrosscorr', false); % Store cross correlation between residual and original matrix?
    storecoefficients = setParam(varargin, 'storecoefficients', false); % Store OMP coefficients?
    
    if MAKEMOVIE; 
        figure; 
        F(1) = struct('cdata',[],'colormap',[]); 
        set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]); 
    end

    r = y;
    Phi = [];
    S = [];
    a = zeros(size(y));
    tvy = tv(y);
    punexp = tv(r)/tvy;
    crosscorr = [];
    while punexp(end) > c && length(S) < maxfeatures;
        
        if MAKEMOVIE
            subplot(3,1,1);
            imagesc(r, [-3 3]); colormap(1 - prgn);
            axis off
            title('Residual', 'FontSize', 20);
            hold on
            for i = 1 : length(S)
                plot([S(i) S(i)], [0 size(y,1)], '-r', 'LineWidth', 1);
            end
            drawnow
            
            
            subplot(3,1,2);
            imagesc(a, [-3 3]); colormap(1 - prgn);
            axis off
            title('Approximation', 'FontSize', 20);
            hold on
            for i = 1 : length(S)
                plot([S(i) S(i)], [0 size(y,1)], '-r', 'LineWidth', 1);
            end
            drawnow
            
            subplot(3,1,3);
            xmax = 25;
            plot(punexp, '-ok', 'LineWidth', 2);
            axis([1 xmax 0 1]);
            hold on
            plot([1 xmax], [c c], '-r');
            title('Proportion of unexplained variance', 'FontSize', 20);
            set(gca, 'FontSize', 18);
            drawnow
            
            F(end+1) = getframe(gcf);
        end
        
        
        % Naively computing the inner product will produce a genes x genes
        % correlation matrix, which may be quite difficult to store in
        % memory.
        if savememory
            [k,z] = max_sum_abs_inner_prod(r',y);
        else
            z = sum(abs(r'*y))/(size(y,2)*size(y,1));
            [~, k] = max(z); 
        end
        
        if storecrosscorr; crosscorr = [crosscorr; z]; end;
        S = [S,k];
        Phi = [Phi, y(:,k)];
        
        b = Phi'*Phi \ Phi'*y;
        a = Phi*b;
        r = y - a;
        
        punexp(end+1) = tv(r)/tvy;
        
        fprintf('Prop. unexplained variance: %0.4f\n', punexp(end));
    end

    stats.S = S;
    if storecoefficients; stats.b = b; else; stats.b = []; end
    punexp(1) = []; % always == 1
    stats.punexp = punexp;
    stats.crosscorr = crosscorr;
    if MAKEMOVIE; F(1) = []; stats.F = F; end

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


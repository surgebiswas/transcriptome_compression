function model = tratrain( y, x, varargin )
% Y = samples x genes
% x = samples x marker_genes 
% 
% Updated tratrain to have lambda multiplied by the size of the data. The
% interpretation of lambda then becomes influence of the prior in terms of
% proportions of available training data. This modification is helpful
% since now the search space for lambda is independent of the size of the
% dataset.

USE1SE = true;

% Standardize variables.
[y,uy,sy] = st(y);
[x,ux,sx] = st(x);



nfold = setParam(varargin, 'nfold', 10);
lambda = setParam(varargin, 'lambda', logspace(-8, 1, 100));

N = size(y,1)*(nfold - 1)/nfold;
I = N*eye(size(x,2)); % Scale by size of training data.

cvind = crossvalind('Kfold',size(y,1),nfold);


mse = zeros(nfold,length(lambda));
xtx = cell(nfold,1);
xty = cell(nfold,1);
for i = 1 : length(lambda)
    fprintf('Iter: %0.0f\tLambda = %0.8f\n', i, lambda(i));
    for j = 1 : nfold
        
        if isempty(xtx{j})
            xtrain = x(cvind ~= j,:);
            ytrain = y(cvind ~= j,:);

            xtest = x(cvind == j, :);
            ytest = y(cvind == j, :);
            
            % Precompute the covariance and cross covariance matrices for
            % each fold.
            xtx{j} = xtrain'*xtrain;
            xty{j} = xtrain'*ytrain;            
        end
        
        bhat = (xtx{j} + lambda(i)*I) \ xty{j}; % Ridge estimate.
        in = isnan(bhat);
        if any(any(in));
            warning('NaN''s found in coefficient estimate. Setting them to 0.');
            bhat(in) = 0;
        end
        
        % Evaluate test error
        yhat = xtest*bhat;
        
        mse(j,i) = mean(mean( (ytest - yhat).^2 ));
    end
end


% Train the full model
[minmse,mind] = min(mean(mse));

% Use lambda value that gives mean MSE that's within 1 standard error of
% the best mean MSE? 
if USE1SE
    se = zeros(1,size(mse,2));
    me = zeros(1,size(mse,2));
    for i = 1 : length(se)
        % remove MSEs from outlying folds before calculating standard
        % errors. This gives a more robust estimate of standard error.
        m = mse(:,i);
        
        sd = std(m);
        m( (m > mean(m) + 2*sd) | (m < mean(m) - 2*sd) ) = [];
        
        se(i) = std(m)/sqrt(length(m));
        me(i) = mean(m);
    end
    mind = find( me < minmse+se, 1,  'last');
end

bstar = (x'*x + lambda(mind)*I) \ x'*y;



% unstandardize back to the original scale.
SY2 = repmat(sy, size(bstar,1), 1);
SX2 = repmat(sx', 1, size(bstar,2));
bstar_us = SY2.*bstar./SX2;


model.mse = mse;
model.lambda_star = lambda(mind);
model.b = bstar_us;
model.b0 = uy - ux*model.b;
model.b_standardized = bstar;

    function [m,um,sm] = st(m)
        % Standardize.
        n = size(m,1);
        
        um = mean(m);
        sm = std(m);
        for k = 1 : size(m,2)
            if sm(k) ~= 0
                m(:,k) = (m(:,k) - um(k))/sm(k);
            else
                m(:,k) = 0;
            end
        end
    end

end


function model = tratrain( y, x, varargin )
% Y = samples x genes
% x = samples x marker_genes 

% Standardize variables.
[y,uy,sy] = st(y);
[x,ux,sx] = st(x);



nfold = setParam(varargin, 'nfold', 10);
lambda = setParam(varargin, 'lambda', logspace(-3, 2, 100));

I = eye(size(x,2));
cvind = crossvalind('Kfold',size(y,1),nfold);

mse = zeros(nfold,length(lambda));
for i = 1 : length(lambda)
    fprintf('Iter: %0.0f\tLambda = %0.6f\n', i, lambda(i));
    for j = 1 : nfold
        
        xtrain = x(cvind ~= j,:);
        ytrain = y(cvind ~= j,:);
        
        xtest = x(cvind == j, :);
        ytest = y(cvind == j, :);
        
        bhat = (xtrain'*xtrain + lambda(i)*I) \ xtrain'*ytrain; % Ridge estimate.
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
[~,mind] = min(mean(mse));
bstar = (x'*x + lambda(mind)*I) \ x'*y;



% unstandardize back to the original scale.
SY2 = repmat(sy, size(bstar,1), 1);
SX2 = repmat(sx', 1, size(bstar,2));
bstar_us = SY2.*bstar./SX2;


model.mse = mse;
model.lambda_star = lambda(mind);
model.b = bstar_us;
model.b0 = uy - ux*model.b;

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


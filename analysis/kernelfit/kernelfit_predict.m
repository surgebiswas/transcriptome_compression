function [ yhat ] = kernelfit_predict( x0, model )

if model.train_standardized
    x0 = standardize(x0, 'mu', model.x_mu, 'std', model.x_sig);
end


yhat = zeros(size(x0,1), size(model.y,2));
for i = 1 : size(x0,1)
    x = model.x;
    y = model.y;
    
    w = model.K(x(:,2:end), x0(i,:), model.lambda);
%     mask = w < 0; %1e-4;
%     
%     x(mask,:) = [];
%     y(mask,:) = [];
%     w(mask) = [];

    xtw = bsxfun(@times, x, w)';

    yhat(i,:) = [1 x0(i,:)]*( xtw*x \ xtw*y );
end

% Return to original scale if requested.
% Use training data mean and standard deviation.
if model.predict_original
    yhat = unstandardize(yhat, model.y_mu, model.y_sig);
end

end


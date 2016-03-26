function [ yhat ] = kernelfit_predict( x0, model )

if model.train_standardized
    x0 = standardize(x0, 'mu', model.x_mu, 'std', model.x_sig);
end


yhat = zeros(size(x0,1), size(model.y,2));
x = model.x;
y = model.y;
for i = 1 : size(x0,1)
    w = model.K(x(:,2:end), x0(i,:), model.lambda);

    if strcmpi(model.method, 'local_regression');
        xtw = bsxfun(@times, x, w)';
        yhat(i,:) = [1 x0(i,:)]*( xtw*x \ xtw*y );
    elseif strcmpi(model.method, 'local_average');
         w = w/sum(w);
         yhat(i,:) = w'*y;
    end
end

% Return to original scale if requested.
% Use training data mean and standard deviation.
if model.predict_original
    yhat = unstandardize(yhat, model.y_mu, model.y_sig);
end

end


function l = kernelfit_loss( y, yhat, model )
% l = kernelfit_loss( y, yhat, model)
% Calculates MSE between y and yhat. 
% 
% y - vector or matrix of target response values.
% yhat - equally sized vector or matrix of predicted response values.
% model - struct output from kernelfit_train. Must contain the fields
%         predict_original. If predict_original is true, then loss is
%         calculated directly between y and yhat. If false, yhat is
%         standardized using mean and standard deviation obtained during
%         training (fields y_mu and y_sig must be specified).


if ~model.predict_original
    % If predictions were not returned to the original scale (i.e. yhat is 
    % standardized), then we must standardize y.
    
    y = standardize(y, 'mu', model.y_mu, 'std', model.y_sig);
end

l = mean(mean((y - yhat).^2));




end


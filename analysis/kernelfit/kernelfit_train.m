function [ model ] = kernelfit_train( x, y, params )
% [ model ] = kernelfit_train( x, y, params )
% "Trains" a local regression model. There is no training really required,
% as the entirety of computation occurs during prediction.
%
% INPUT
% x: An [n-observations x p-features] design matrix.
% y: An [n-observations x r-responses] response matrix.
% params: A cell array vector of length 2. 
%   - params{1} should contain a kernel function handle. The kernel function
%       should have three arguments.
%           x: An [n-observations x p-features] training design matrix. 
%           x0: A single [1 x p-features] query point.
%           lambda: The kernel width parameter.   
%       As an example, the tricube kernel function can be defined as:
%       K = @(x,x0,lambda) (1 - min(abs(bsxfun(@minus,x,x0)/lambda), 1).^3).^3 
%
%   - params{2} should contain lambda, the kernel width (see above).
%
%   - params{3} Train using a standardized model? [true | false]  
%   
%   - params{4} Predict in the original scale? [true | false]  
% 
% Note: Written to be used with cvalidate.m and cvalidate_tune.m


model.K = params{1}; % Kernel
model.lambda = params{2}; % Kernel width. 
model.train_standardized = params{3}; % train on a standardized model?
model.predict_original = params{4}; % make predictions in original scale? 

if model.train_standardized
    [model.x, model.x_mu, model.x_sig] = standardize(x); 
    model.x = [ones(size(x,1),1), model.x];
    [model.y, model.y_mu, model.y_sig] = standardize(y);
else
    model.x = [ones(size(x,1),1), x];
    model.y = y;
end



end


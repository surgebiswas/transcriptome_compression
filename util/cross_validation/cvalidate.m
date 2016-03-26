function [ cvperf_mean, cvperf_se, cvperf ] = cvalidate( x, y, trainfun, predfun, varargin  )
% [ cvperf_mean, cvperf_se, cvperf ] = cvalidate( x, y, trainfun, predfun, varargin  )
% cvalidate - General interface to perform cross validation.
% Surojit Biswas
%
% INPUT
% x: An [n-observations x p-features] design matrix.
%
% y: An [n-observations x r-responses] target matrix.
% 
% trainfun: A function handle to a function with three required inputs:
%   - xtrain: An [n-observations x p-features] training design matrix.
%   - ytrain: An [n-observations x r-responses] training target matrix.
%   - params: A variable containing parameters required for training. If
%       not required, this argument should be left empty ([]).
%   
%   trainfun should output a single variable which contains model
%   parameters that can be passed to predfun (see predfun argument below).
% 
% predfun: A function handle to a function with two inputs:
%   - xtest: An [n-observations x p-features] test design matrix.
%   - model: A model structure output by trainfun that contains the
%       required model parameters for testfun to make a prediction given 
%       xtest.
%
% OPTIONAL INPUTS
% crossvalind: A [size(x,1) x 1] vector of cross validation indices.
%   Default = crossvalind('Kfold', size(x,1), 10)
% 
% training_params: A variable (e.g. struct, cell array, vector) containing 
%   training parameters required for training by trainfun. Leave empty if
%   irrelevant.
%   Default = []
% 
% loss_fun: A function handle for a loss function.
%   Default = @(ytest,yhat) mean(mean((ytest - yhat).^2)) (mean squared error)
%
%
% OUTPUT
% cvperf_mean: Mean performance (loss) across all folds. 
% cvperf_se: Standard error of the performance across all folds.
% cvperf: Cross validation performance for each fold. 

DEBUG_VERBOSE = true;

cvi = setParam(varargin, 'crossvalind', crossvalind('Kfold', size(x,1), 10));
params = setParam(varargin, 'training_params', []);
lossfun = setParam(varargin, 'loss_fun', @(ytest,yhat,s) mean(mean((ytest - yhat).^2)));


cvperf = zeros(max(cvi),1);
for i = 1 : length(cvperf)
    xtrain = x(cvi~=i,:);
    ytrain = y(cvi~=i,:);

    xtest = x(cvi==i,:);
    ytest = y(cvi==i,:);
    
    model = trainfun(xtrain,ytrain,params);
    yhat = predfun(xtest,model);
    
    cvperf(i) = lossfun(ytest,yhat);
    
    if DEBUG_VERBOSE
        fprintf('[%s] CV index: %0.0f\n', mfilename, i);
    end
end

cvperf_mean = mean(cvperf);
cvperf_se = std(cvperf)/sqrt(length(cvperf));


end


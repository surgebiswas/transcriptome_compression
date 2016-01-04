function cvstats = cvalidate_tune( x, y, trainfun, predfun, params, varargin )
% cvstats = cvalidate_tune( x, y, trainfun, predfun, params, varargin )
%
% cvalidate_tune - Uses cross-validation test-set accuracy to tune
% hyperparameters associated with a supervised learning algorithm.
%
% INPUT
% x: An [n-observations x p-features] design matrix.
%
% y: An [n-observations x r-responses] target matrix.
%
% trainfun: A function handle to a function with three required inputs:
%   - xtrain: An [n-observations x p-features] training design matrix.
%   - ytrain: An [n-observations x r-responses] training target matrix.
%   - param: A [1 x #-parameters] cell array containing hyperparameters 
%       required for training. Note that this format requirement is more
%       strict than it is for cvalidate.m. For example, suppose trainfun
%       requires a mean vector and covariance matrix as hyperparameters.
%       Then trainfun may be called as trainfun(xtrain,ytrain, {mu,
%       covar}), where mu is a d-dimensional vector and covar is a dxd
%       covariance matrix.
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
% params: A [1 x #-hyperparameters] cell aray. Each cell within this cell
%   array should contain either a numeric vector or another cell array
%   vector defining a schedule of parameter values to perform the cross
%   validation search over. For example, if params{2} = logspace(-6,3,100),
%   then this would define a logrithmically spaced schedule of length 100 
%   (from 10^-6 to 10^3) for the second hyperparameter. 
%
% OPTIONAL INPUTS
% crossvalind: A [size(x,1) x 1] vector of cross validation indices.
%   Default = crossvalind('Kfold', size(x,1), 10)
% 
% loss_fun: A function handle for a loss function.
%   Default = @(ytest,yhat) mean(mean((ytest - yhat).^2)) (mean squared error)
%
% OUTPUT
% cvstats: A structure with the following fields:
%   - param_idx = A [q x #-hyperparameters] matrix of indices used to
%   define hyperparameter combinations used during training. Here q is
%   equal to prod(cellfun(@length,params)).
%   For example, if param_idx(4,:) = [1 2 2] then during the 4th iteration
%   of the cross validation search, training occurred with [params{1}(1),
%   params{2}(2), params{3}(2)].
%
%   - cvperf_mean = [q x 1] vector of average test set loss for each
%   hyperparameter combination explored.
% 
%   - cvperf_se = [q x 1] vector of the test set loss standard error for each
%   hyperparameter combination explored.
%
%   - param_idx_best_perf = [1 x #-hyperparameters] index vector specifying
%   the combination of hyperparameters that give the best test set
%   performance. This is a row of param_idx.
%
%   - param_best_perf = [1 x #-hyperparameters] cell array containing the
%   hyperparameter combination that gives the best test set performance.
%
%   - full_model = trainfun(x,y,cvstats.param_best_perf), where x and y 
%   represent the full training data. 

cvi = setParam(varargin, 'crossvalind', crossvalind('Kfold', size(x,1), 10));
lossfun = setParam(varargin, 'loss_fun', @(ytest,yhat) mean(mean((ytest - yhat).^2)));

% Generate evaluation indices over parameters.
pidx = zeros(prod(cellfun(@length,params)), length(params));
pidx(1,:) = ones(1, size(params,2));
for j = 2 : size(pidx,1)
    pidx(j,:) = get_next_param_idx(params, pidx(j-1,:));
end

% Perform cross validation.
cvperf_mean = zeros(size(pidx,1),1);
cvperf_se = zeros(size(pidx,1),1);
for j = 1 : size(pidx,1)
    [cvperf_mean(j), cvperf_se(j)] = cvalidate(x, y, trainfun, predfun, ...
        'training_params', build_parameter_vector(params, pidx(j,:)), 'crossvalind', cvi, ...
        'loss_fun', lossfun);
end

% Select the best parameter value.
[~,mind] = min(cvperf_mean); % assumes performance is loss.

cvstats.param_idx = pidx;
cvstats.cvperf_mean = cvperf_mean;
cvstats.cvperf_se = cvperf_se;
cvstats.param_idx_best_perf = pidx(mind,:);
cvstats.param_best_perf = build_parameter_vector(params, cvstats.param_idx_best_perf);
cvstats.full_model = trainfun(x,y,cvstats.param_best_perf);



    function pv = build_parameter_vector(params, pvi)
        pv = cell(1,length(pvi));
        for i = 1 : length(pvi)
            if isnumeric(params{i}) || islogical(params{i})
                pv{i} = params{i}(pvi(i));
            else
                pv{i} = params{i}{pvi(i)};
            end
        end
    end

    function current_idx = get_next_param_idx(params, current_idx)
        for i = length(current_idx) : -1 : 1
            if current_idx(i) + 1 <= length(params{i})
                current_idx(i) = current_idx(i) + 1;
                break
            else
                current_idx(i) = 1;
            end
        end
    end


end


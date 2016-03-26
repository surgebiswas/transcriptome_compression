function [ model ] = ridgefit_train( x, y, params )
% [ model ] = ridgefit_train( x, y, params )
%
% Trains a ridge regression model. 
%
% INPUT
% x: An [n-observations x p-features] design matrix.
% y: An [n-observations x r-responses] response matrix.
%
% - params{1} = L2 penalty (lambda). Note that during pseudoinverse
%       calculation, lambda is scaled by number of observations (i.e.
%       lambda_actual = size(x,1)*params{1}. 
% - params{2} = Compute x'*x and x'*y or use precomputed versions? Using
%       precomputed versions is more efficient for cross validation
%       purposes. For training a model de novo, set params{2} to false.
% - params{3} = If params{2} is true, params{3} needs to be a
%       ridgefit_precompute object with fields populated as follows:
%       - precompute: An [f x 2] cell array contain precomputed versions of
%       x'*x and x'*y in columns 1 and 2 respectively.
%       - itr: Iterator index that specifies which row of precompute to
%       use. itr is incremented after every call of ridgefit_train and is
%       reset to 1 after it equals size(precompute,1).
%
%       Using params{3} is useful for cross validation based tuning of the
%       L2 penalty. The ridge estimate of the regression coefficients is
%       given by (x'*x + lambda*I) / x'*y. If cross validation indices are
%       set before time x'*x and x'*y can be computed for each training
%       fold separately one time. These results can be stored in params{3}
%       as a num_folds x 2 cell array for reuse with every lambda examined.
%
%       IMPORTANT: When precomputing, make sure a column of ones is added 
%       for the intercept term.
% 
% Note: written to be used with cvalidate.m and cvalidate_tune.m

lam = size(x,1)*params{1}*eye(size(x,2) + 1); % the plus 1 is for an intercept column.
if params{2}  
    
    p = params{3};
    xtx = p.precompute{p.itr,1};
    xty = p.precompute{p.itr,2};
    
    % If x'*x was computed with a column of ones, then x'*x(1,1) should
    % equal size(x,1), the number of training observations.
    if (xtx(1,1) == size(x,1))
       
%         assert(xtx(1,1) == size(x,1), ['Precomputed x''*x is not consistent ', ...
%             'with the supplied design matrix. Either an intercept column of ', ...
%             'ones was not added during precomputation or the iterator index ', ...
%             'has been offset.']);

        model.b = (xtx + lam) \ xty;
        model.precompute_used = true;

        % Update the iterator.
        p.itr = mod(p.itr+1,size(p.precompute,1));
        if p.itr == 0; p.itr = size(p.precompute,1); end
    else
        model.b = estimate_denovo([ones(size(x,1),1), x],y,lam);
        model.precompute_used = false;
    end
else
    model.b = estimate_denovo([ones(size(x,1),1), x],y,lam);
    model.precompute_used = false;
end


    function b = estimate_denovo(x,y,lam)
        b = (x'*x + lam) \ x'*y;
    end

end


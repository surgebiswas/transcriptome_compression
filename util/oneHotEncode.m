function [ x, fnames ] = oneHotEncode( xd, varargin )
% Converts design matrix dataset xd into a numerical design matrix directly
% useable for supervised learning applications.
%
% xd =  [n x p] design matrix dataset. Can contain numerical and categorical
%       predictors. Categorical predictors are assumed to be represented as
%       strings.

rmli = setParam(varargin, 'rmli', true);
addIntercept = setParam(varargin, 'intercept', true);



% Predictor names
vn = get(xd, 'VarNames');

% Determine categorical predictors
iscateg = false(1,size(xd,2));
for i = 1 : size(xd,2)
    iscateg(i) = iscell( eval(['xd.', vn{i}, '(1)']) );
end

x = [];
fnames = [];
for i = 1 : size(xd,2)
    if iscateg(i)
        [xs, fs] = buildSubMatrix(eval(['xd.', vn{i}]), vn{i});
    else
        xs = eval(['xd.', vn{i}]);
        fs = vn(i);
    end
    
    x = [x,xs];
    fnames = [fnames,fs];
end

if addIntercept
    x = [ones(size(x,1),1),x];
    fnames = [{'Intercept'}, fnames];
end



function [x, f] = buildSubMatrix(column, varname)
        o = double(rmli);
        u = flipud(unique(column));
        
        x = zeros(length(column), length(u) - o);
        
        f = u(1:end-o)';
        for q = 1 : length(u) - o
            x(:,q) = double(strcmpi(column, u{q}));
            f{q} = [varname, '_', f{q}];
        end
    end



end


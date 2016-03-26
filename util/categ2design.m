function [ X, features ] = categ2design( cat, addIntercept, varargin )
% categ2design - Creates a numerical design matrix representation of a
% categorical matrix. 
% 
% INPUT:
% cat           =   [n x p] cell array or dataset containing categories for each of n
%                   samplse
%
% OUTPUT:
% X             =   Numerical design matrix representation of cat.
% features      =   Column labels for X which describes what each predictor
%                   is.


if strcmpi(class(cat), 'dataset')
    vn = get(cat, 'VarNames');
    
    if ~isempty(get(cat, 'ObsNames'))
        cat = dataset2cell(cat);
        cat = cat(2:end,2:end); % remove row and column names.
        
    else
        cat = dataset2cell(cat);
        cat = cat(2:end,:);
    end
end

rmli = getVarargin(varargin, 'rmli'); % remove extra covariate vector that reduces rank?
if isempty(rmli)
    rmli = true;
end

s = getVarargin(varargin, 'sparse');
if isempty(s)
    s = false;
end

trace = getVarargin(varargin, 'trace');
if isempty(trace)
    trace = false;
end

rowreduce = getVarargin(varargin, 'rref');
if isempty(rowreduce)
    rowreduce = false;
end

interaction = getVarargin(varargin, 'interaction');
if isempty(interaction)
    interaction = false;
end


if iscell(cat)
    treatAsString = ischar(cat{1,1});
else
    treatAsString = false;
end

n = size(cat,1);
p = size(cat,2);

X = [];
features = [];
if interaction; blockSize = zeros(1,p); end
for k = 1 : p
    if trace; fprintf('One hot encoding column %0.0f/%0.0f ... ', k, p); end
    
    if ischar(cat{1,k})
        [x, f] = buildSubMatrix(cat(:,k), s);
    else
        % assume its numerical
        x = cell2mat(cat(:,k));
        f = vn(k);
    end
    X = [X, x];
    features = [features, f];
    blockSize(k) = length(f);
    if trace; fprintf('Done.\n'); end
end

% Build interaction terms if requested.
if interaction
    totalInteractions = length(blockSize)*(length(blockSize) - 1)/2;
    blocksBuiltIdx = 0;
    for k = 1 : length(blockSize)
        for l = k + 1 : length(blockSize)
            if trace; fprintf('One hot encoding interactions between column %0.0f (# unique values = %0.0f) and column %0.0f (# unique values = %0.0f) ... ', k, blockSize(k), l, blockSize(l)); end
            abk = sum(blockSize(1:k-1));
            interactoridx =  (abk + 1) : (abk + blockSize(k));

            abl = sum(blockSize(1:l-1));
            interacteeidx = (abl + 1) : (abl + blockSize(l));

            % Allocate space 
            if s
                x = spalloc(size(X,1), length(interactoridx)*length(interacteeidx), nnz(X(:, interactoridx)) + nnz(X(:,interacteeidx)));
            else
                x = zeros(size(X,1), length(interactoridx)*length(interacteeidx));
            end

            idx = 1;
            for ii = 1 : length(interactoridx)
                fprintf('Encoding sub-block %0.0f.\n', ii); 
                for jj = 1 : length(interacteeidx)
                    x(:,idx) = X(:,interactoridx(ii)).*X(:,interacteeidx(jj));
                    idx = idx + 1;
                end
                
            end

            X = [X, x];
            blocksBuiltIdx = blocksBuiltIdx + 1;
            if trace; fprintf('Done. %0.0f/%0.0f interactions encoded.\n', blocksBuiltIdx, totalInteractions); end
        end
    end
end

% Remove all zero columns
X(:, sum(X) == 0) = [];


if addIntercept
    X = [ones(n,1), X];
    if treatAsString
        features = [{'Intercept'}, features];
    else
        features = [1 features];
    end
end
%assert(size(X,2) == length(features));

% Remove low rank columns
if rowreduce
[~, licols] = rref(X);
X = X(:, licols);
features = features(licols);
end

%assert(size(X,2) == length(features));



    function [x, f] = buildSubMatrix(column, buildSparse)
        if rmli
            o = 1;
        else
            o = 0;
        end
        u = unique(column);
        
        if buildSparse
            x = spalloc(length(column), length(u) - o, length(column));
        else
            x = zeros(length(column), length(u) - o);
        end
        
        f = u(1:end-o)';
        for i = 1 : length(u) - o
            if treatAsString
                x(:,i) = double(strcmpi(column, u{i}));
            else
                x(:,i) = double(column == u(i)); % Assume this is numerical matrix.
            end
        end
    end


end


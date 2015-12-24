function [ cvi ] = kfoldcrossvalindbygroup( k, g )
% Constructs kfold cross validation indices; however, samples with the same 
% grouping variable are assigned to the same fold.
% 
% k = number of desired folds.
% g = vector or cell array of grouping IDs. g(i) is a number or g{i} is a 
%   string that specifies the group observation i is in. length(g) = number 
%   of observations.

ug = unique(g);

c = crossvalind('Kfold', length(ug), k);

cvi = zeros(length(g),1);
for i = 1 : length(ug)
    
    if iscell(ug)
        cvi( strcmp(g, ug{i}) ) = c(i);
    else % assume its double. Change later if needed.
        cvi( g == ug(i) ) = c(i);
    end
end




end


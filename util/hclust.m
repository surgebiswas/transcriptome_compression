function [ ri, ci ] = hclust( x )
% Hierarchical agglomerative clustering. Functionally equivalent to
% clustergram(x, 'OptimalLeafOrder', false);
%
% x = an n x p matrix. e.g. a samples x genes matrix of standardized
% expression values.
%
% OUTPUT:
% ri = permutation vector for row indices.
% ci = permutation vector for column indices.


ri = dendroperm(linkage(x, 'average'), 0);
ci = dendroperm(linkage(x', 'average'), 0);




end


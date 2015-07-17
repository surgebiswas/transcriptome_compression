function [ ri, ci ] = hclust( x, varargin)
% Hierarchical agglomerative clustering. Functionally equivalent to
% clustergram(x, 'OptimalLeafOrder', false);
%
% x = an n x p matrix. e.g. a samples x genes matrix of standardized
% expression values.
%
% OUTPUT:
% ri = permutation vector for row indices.
% ci = permutation vector for column indices.

useksvd = setParam(varargin, 'useksvd', []);

if isempty(useksvd)
    fprintf('Clustering rows ... ');
    ri = dendroperm(linkage(x, 'average'), 0);
    fprintf('Done.\n');
    
    fprintf('Clustering columns ... ');
    ci = dendroperm(linkage(x', 'average'), 0);
    fprintf('Done.\n');
else
    fprintf('Computing SVD (k = %0.0f) ... ', useksvd);
    [U,S,V] = svdsecon(x, useksvd);
    fprintf('Done.\n');
    
    fprintf('Clustering rows ... ');
    ri = dendroperm(linkage(U*S, 'average'), 0);
    fprintf('Done.\n');
    
    fprintf('Clustering columns ... ');
    ci = dendroperm(linkage(V*S, 'average'), 0); % Note V*S = (S*V')'
    fprintf('Done.\n');
end
    
    






end


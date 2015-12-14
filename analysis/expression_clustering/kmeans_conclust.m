function [ c, csoft, cmat ] = kmeans_conclust( x, varargin )

KMAX = 20;
NCLUSTEVALREPS = 5;
nrep = setParam(varargin, 'nrep', 100);

optKs = [];
for i = 1 : NCLUSTEVALREPS
    eva1 = evalclusters(x, 'kmeans', 'CalinskiHarabasz', 'klist', 1:KMAX);
    eva2 = evalclusters(x, 'kmeans', 'DaviesBouldin', 'klist', 1:KMAX);
    optKs = [optKs, eva1.OptimalK];
    optKs = [optKs, eva2.OptimalK];
end



optK = floor(mean(optKs));

cmat = zeros(size(x,1), nrep);
for i = 1 : nrep
    cmat(:,i) = kmeans(x,optK, 'start', 'uniform');
end

[c, csoft] = conclust(cmat);

% Note that size(csoft,2) and max(c) may be different values after removing
% zero cluster labels. 
c = remove_zero_cluster_labels(c);


end


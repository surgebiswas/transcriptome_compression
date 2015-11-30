function [ c, csoft, cmat ] = kmeans_conclust( x, varargin )

KMAX = 20;
nrep = setParam(varargin, 'nrep', 100);

eva1 = evalclusters(x, 'kmeans', 'CalinskiHarabasz', 'klist', 1:KMAX);
eva2 = evalclusters(x, 'kmeans', 'DaviesBouldin', 'klist', 1:KMAX);

optK = floor(mean([eva1.OptimalK, eva2.OptimalK]));

cmat = zeros(size(x,1), nrep);
for i = 1 : nrep
    cmat(:,i) = kmeans(x,optK, 'start', 'uniform');
end

[c, csoft] = conclust(cmat);




end


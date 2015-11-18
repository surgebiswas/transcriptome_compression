function Ensemble_labels = conclust( c )
% Consensus clustering
%
% Multiple cluster index assignments. n-observations x k-clusterings matrix
% of cluster indices. See apcluster_boot

PRAMALAP = 0.000001;

nclust = round(mean(max(c,[],1))) - 3;

palpha = ones(nclust,1)/nclust;
pbeta = cell(1,size(c,2));
for i = 1 : length(pbeta)
    pbeta{i} = ones(nclust,max(c(:,i)))/max(c(:,i));
end

[phiAll,gammaAll,resultAlpha,resultBeta] = learnBCE(c, palpha, ...
    pbeta, PRAMALAP, max(c,[],1));

% Obtain the cluster assignments from BCE
Ensemble_labels=zeros(size(c,1),1);
for index=1:size(c,1)
    wtheta(:,index)=gammaAll(:,index);
    bb=find(wtheta(:,index)==max(wtheta(:,index)));
    Ensemble_labels(index)=bb(1);
end



end


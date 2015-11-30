function [ c ] = apcluster_boot( x, b )
% Bootstrap apcluster
% 
% x = n-observations x p-features data matrix
% b = number of bootstraps to perform.
%
% c = n x 1 vector of cluster assignments

K = round(0.01*size(x,1));

pd = squareform(pdist(x));
c = [];
for i = 1 : b
    fprintf('Running affinity propagation ... Bootstrap idx = %0.0f\n', i);
    bind = randsample(size(x,1),size(x,1), true);
    
    pdb = pd(bind,bind);
    
    idx = apcluster(pdb,1);
    
    % Replace cluster indices with starting positive integers.
    ui = unique(idx);
    for j = 1 : length(ui)
        idx(idx == ui(j)) = j;
    end
    
    % We now have cluster indices for all observations that made it into
    % the bootstrap. We now to need to get cluster indices over the whole
    % dataset.
    % 
    % Impute these indices using knn searches.
    idxfull = nan(size(idx));
    idxfull(bind) = idx;
    
    ubind = unique(bind);
    ubind_idx = idxfull(ubind);
    for j = 1 : length(idxfull)
        if isnan(idxfull(j))
            knnidx = knnsearch(x(ubind,:), x(j,:), 'K', K);
            
            ucid = unique(ubind_idx(knnidx));
            if length(ucid) > 1
                ccounts = hist(ubind_idx(knnidx), ucid);

                [~,mind] = max(ccounts);
                if mind < 1 || mind > length(ucid)
                    keyboard;
                end
                imp_idx = ucid(mind);
                idxfull(j) = imp_idx;
            else
                idxfull(j) = ucid;
            end
        end
    end
    
    c = [c, idxfull];
end




end


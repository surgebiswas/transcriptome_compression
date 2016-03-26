function tradict_train( lY, trainfun, predfun, paramsetupfun, varargin )
%UNTITLED18 Summary of this function goes here
%   Detailed explanation goes here

control.rep_target_idx = setParam(varargin, 'rep_target_idx', 'random');
D = setParam(varargin, 'distance_matrix', []); % Must be in NON-squareform M*(M-1)/2 long vector.

% Initialize
sY = standardize(lY);
progress.cidx = ones(size(lY,2),1);
progress.cidx_star = 1;

if isempty(D)
    progress.D = pdist(sY');
else
    progress.D = D;
end

while not_converged
    progress.cavg = compute_cluster_averages(sY, progress);
    idx = find_best_gene_idx(sY, progress);
end

    function idx = find_best_gene_idx(lY, progress)
        
        
        
        for i = 1 : size(cavg,2)
            x = [progress.selected_features, cavg(:,i)];
            
            
        end
        
        
    end




end


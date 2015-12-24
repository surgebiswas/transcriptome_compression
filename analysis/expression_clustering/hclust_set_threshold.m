function thresh =  hclust_set_threshold( x, varargin )
% Uses 
% Assumes clustering is on the rows of x.
% Cross validation indices define a partition over the columns of x.
%
% e.g. x = [genes x samples]. 

z = setParam(varargin, 'linkage_output', []);
cvi = setParam(varargin, 'crossvalind', crossvalind('Kfold', size(x,2), 10));

assert(~isempty(z), 'Currently require user supplied output from linkage.m');

cutoffs = sort(unique(z(:,end))) - 1e-6; cutoffs(1) = 0;
testperformance = nan(size(cutoffs));


% Use bisection (binary search) to explore cutoff space.
not_converged = true;
bisect_left = 1;
bisect_right = length(cutoffs);
while not_converged
    
end


    function evaluate_test_set_performance(x, z, cvi, c)
        cidx = cluster(z, 'cutoff', c, 'criterion', 'distance');
        if max(cidx) > 
        
        for i = 1 : max(cvi)
            
        end
    end


end


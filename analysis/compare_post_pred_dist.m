function [ ci_sched, act_contained ] = compare_post_pred_dist( ss, post_smpl )

% ss = samples x features matrix of actual values.
% post_smpl = samples x 1 cell array of draws form the posterior

ci_sched = 0:0.05:1;
act_contained = zeros(1, length(ci_sched));
for k = 1 : length(ci_sched)
    disp(k);
    L = zeros(size(ss));
    U = L;
    cred_interval_size = ci_sched(k);
    l = (1 - cred_interval_size)/2; r = 1 - l;
    for i = 1 : size(ss,1)

        p = prctile(post_smpl{i}, 100*[l, r]);
        L(i,:) = p(1,:);
        U(i,:) = p(2,:);
    end
    
    act_contained(k) = sum(ss(:) < U(:) & ss(:) > L(:))/numel(ss);
end

end


function c = robust_cluster_relabel( x, c )
% x = n-observations x p-features dataset that's been clustered.
% c = n-observations vector of cluster labels. Cluster labels are assumed
%       to take values from 1,2,...,k

RDCUTOFF = 3; % In number of standard deviations.
MAXITR = 20;

% Remove cluster indices with zero representation. Just relabel higher
% index clusters to take the empty lower index ones.
c = remove_zero_cluster_labels(c);


C = max(c);
cold = zeros(size(c,1),1);
means = cell(max(c),1);
covs = cell(max(c),1);
rd = zeros(size(c,1),1);
itr = 0;
while any(cold ~= c) && itr < MAXITR
    % Store old label vector for comparison.
    cold = c;
    
    % Compute robust estimates of the mean and covariance for each cluster.
    for i = 1 : C
        rew = mcdcov(x(c == i,:), 'plots', 0);
        means{i} = rew.center;
        covs{i} = rew.cov;
        rd(c == i) = rew.rd;
    end
    
    % Take points with large robust distances from their cluster centers,
    % evaluate their likelihood wrt to other clusters, and reassign if
    % another cluster gives a higher likelihood for that point.
    for i = 1 : size(c,1)
        if rd(i) > RDCUTOFF
            
            % Calculate log-likelihoods.
            loglike = zeros(1,C);
            for j = 1 : C
                loglike(j) = logmvnpdf(x(i,:), means{j}, covs{j});
            end
            
            [~,c(i)] = max(loglike);
        end
    end
    
    fprintf('Percent of labels re-assigned: %0.2f\n', 100*sum(cold ~= c)/length(c));
    itr = itr + 1;
    
    if true
        clf
        hold on
        rng('default');
        for i = 1 : max(c);
            plot3(x(c == i,1), x(c == i,2), x(c == i,3), '.', 'Color', rand(1,3));
        end
        drawnow;
    end
end




end


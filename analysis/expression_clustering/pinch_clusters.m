function [ c ] = pinch_clusters( x, c, nmin )
% Pinches clusters with few datapoints into larger clusters.
% 
% x = n-observations x p-features dataset that's been clustered.
% c = n-observations vector of cluster labels. Cluster labels are assumed
%       to take values from 1,2,...,k
% nmin = minimum number of data points allowed per cluster. 

t = tabulate(c);
while any(t(:,2) < nmin)
    [~,pci] = min(t(:,2)); % cluster index with the fewest data points.
    
    % Compute robust estimates of the mean and covariance for each cluster.
    C = max(c);
    means = cell(C,1);
    covs = cell(C,1);
    for i = 1 : C
        rew = mcdcov(x(c == i,:), 'plots', 0);
        means{i} = rew.center;
        covs{i} = rew.cov;
    end
    
    % Re-assign each point in the 'pci' cluster.
    pcimask = c == pci;
    for i = 1 : length(pcimask)
        if pcimask(i)
            % Calculate log-likelihoods.
            loglike = zeros(1,C);
            for j = 1 : C
                loglike(j) = logmvnpdf(x(i,:), means{j}, covs{j});
            end
            
            loglike(pci) = -Inf;
            [~,c(i)] = max(loglike);
        end
    end
    
    % Re-assign labels for higher number clusters.
    c(c > pci) = c(c > pci) - 1;
    t = tabulate(c);
end


end


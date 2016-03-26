function [ cidx ] = cluster_by_wcvar( x, z, wcvar_cutoff, varargin )
% cidx = cluster_by_wcvar( z, wcvar_cutoff )
%
% Determines cutoff height for a hierarchically clustered dataset by
% setting it to the maximum height such that the within cluster variation
% (wcvar) is less than wcvar_cutoff.
%
% x = observations x features data matrix. Assumes clustering was done on
% observations (rows).
% z = dendrogram output from linkage.m
% wcvar_cutoff = cutoff for within cluster variance. Specified as a
% proportion of total variance (e.g. 0.1). 
%
% cidx = vector of cluster indices.


converged_var_gap = setParam(varargin, 'converged_var_gap', 0.02);

totalvar = within_cluster_variance(x,ones(size(x,1),1));

% Cutoff threshold is found using sequential linear interpolation.
% Let yc = wcvar_cutoff, y1 = wcvarprev, y2 = wcvar (current), x1 =
% numcutprev, x2 = numcut (current). 
%
% Eqn of a line in point slope form:
% y - y2 = [(y2-y1)/(x2-x1)]*(x - x2) 
% 
% Want to use this equation to estimate for what x, y = yc. Call this point
% x_new
% So, we have
% yc - y2 = [(y2-y1)/(x2-x1)]*(x_new - x2) 
%
% Solving ...
% [(yc-y2)(x2-x1)/(y2-y1)] + x2 = x_new.

numcut = floor(prod(log(size(x))));
numcutprev = 0;
wcvarprev = totalvar;
wcvar_cutoff = wcvar_cutoff*totalvar; % Update the cutoff to be an absolute value.
wcvar = 0;
while abs(wcvar - wcvar_cutoff) > converged_var_gap
    cidx = cluster(z, 'maxclust', numcut);
    wcvar = within_cluster_variance(x,cidx);
    
    if true; 
        fprintf('WC-Var/Total-Var: %0.4f\tNum. clusters: %0.0f\n', wcvar/totalvar, numcut);
    end
    
    % Update cluster number.
    numcut_new = round((wcvar_cutoff - wcvar)*(numcut - numcutprev)/(wcvar - wcvarprev) + numcut);
    
    % Store previous cluster number and wc variance.
    numcutprev = numcut;
    wcvarprev = wcvar;
    
    % Ensure proper boundaries.
    numcut = min(size(x,1), max(1, numcut_new)); 
    
    
end


    function wcvar = within_cluster_variance(x,cidx)
        % cidx is a cluster index vector defining a clustering over the
        % rows of x.
        %
        % Assumes unique(cidx) = [1, 2, ... max(cidx)]
        
        cvars = zeros(max(cidx),1);
        for i = 1 : length(cvars)
            cvars(i) = mean(var(x(cidx == i,:),1,1), 2); % compute variance with normalization N and along dimension 1.
        end
        wcvar = mean(cvars);
    end

end


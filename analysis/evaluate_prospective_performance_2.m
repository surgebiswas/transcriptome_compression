function perfstats = evaluate_prospective_performance_2( lY, x, qt, predfun, trainfun, params )

% leave one submission out CV
usubs = unique(qt.Submission);
cvi = kfoldcrossvalindbygroup(length(usubs), qt.Submission);
lY_hat = zeros(size(lY));
for i = 1 : max(cvi)
    fprintf('fold: %0.0f\n', i);
    xtrain = x(cvi ~= i, :);
    ytrain = lY(cvi ~= i, :);
    
    model = trainfun(xtrain, ytrain,params);
    
    xtest = x(cvi == i, :);
    lY_hat(cvi == i, :) = predfun(xtest, model);
end
    
    [lY_adj, sn] = subadjust(lY, qt.Submission);
    [lY_hat_adj] = subadjust(lY_hat, qt.Submission);


    %perfstats.density_plot_global = pred_v_actual_density_plot(standardize(lY), standardize(lY_hat));
    % Don't consider single samples. 
    perfstats.density_plot_subadj = pred_v_actual_density_plot(lY_adj(sn >= 6,:), lY_hat_adj(sn >= 6,:));
    perfstats.lY = lY;
    perfstats.lY_hat = lY_hat;
    
    
    
end


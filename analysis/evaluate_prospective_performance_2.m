function perfstats = evaluate_prospective_performance_2( lY, x, qt, predfun, trainfun, params )

% leave one submission out CV
usubs = unique(qt.Submission);
cvi = kfoldcrossvalindbygroup(length(usubs), qt.Submission);
lY_hat = zeros(size(lY));
for i = 1 : max(cvi)
    xtrain = x(cvi ~= i, :);
    ytrain = lY(cvi ~= i, :);
    
    model = trainfun(xtrain, ytrain,params);
    
    xtest = x(cvi == i, :);
    lY_hat(cvi == i, :) = predfun(xtest, model);
end
    
    lY_adj = subadjust(lY, qt.Submission);
    lY_hat_adj = subadjust(lY_hat, qt.Submission);


    perfstats.density_plot_global = pred_v_actual_density_plot(lY, lY_hat);
    perfstats.density_plot_subadj = pred_v_actual_density_plot(lY_adj, lY_hat_adj);
    
    
    
end


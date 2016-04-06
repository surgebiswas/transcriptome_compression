function results = evaluate_prospective_performance_2( lY, tids, sets, qt, predfun, trainfun, params )

% leave one submission out CV
rng('default');
progfile = 'evaluate_prospective_performance_2_progress_tmp.mat';
if exist(progfile, 'file')
    load(progfile)
    istart = idx + 1;
else
    istart = 1;
    usubs = unique(qt.Submission);
    cvi = kfoldcrossvalindbygroup(length(usubs), qt.Submission);
    target_proc = zeros(size(lY,1), length(sets));
    pred_proc = zeros(size(lY,1), length(sets));
end

lf = fopen('log.txt', 'w+');
for i = istart : max(cvi)
    fprintf('Fold: %0.0f\n', i);
    fprintf(lf, 'Fold: %0.0f\n', i);
    idx = i; %track the loop variable for progress saving.
    
    % Training
    model = tradict_train( lY(cvi ~= i,:), tids, sets, trainfun, params );
    
    % Build the gene and gene-set target.
    sY_target = standardize(lY(cvi == i,:), 'mu', model.train_mu, 'std', model.train_sig);
    target_proc(cvi == i, :) = sY_target*model.geneset.coef;
    
    % Formulate the prediction
    pred_proc(cvi == i, :) = predfun(lY(cvi == i,model.S), model.fit);
    
    if mod(i,5) == 0
        save(progfile);
    end
end
fclose(lf);
results.cvi = cvi;
results.target_proc = target_proc;
results.pred_proc = pred_proc;

delete(progfile);
    
%     [lY_adj, sn] = subadjust(lY, qt.Submission);
%     [lY_hat_adj] = subadjust(lY_hat, qt.Submission);
% 
% 
%     %perfstats.density_plot_global = pred_v_actual_density_plot(standardize(lY), standardize(lY_hat));
%     % Don't consider single samples. 
%     perfstats.density_plot_subadj = pred_v_actual_density_plot(lY_adj(sn >= 6,:), lY_hat_adj(sn >= 6,:));
%     perfstats.lY = lY;
%     perfstats.lY_hat = lY_hat;
    
    
    
end


function results = evaluate_prospective_performance_3( Y, o, tids, sets, qt, nfolds )
% Uses Poisson MVN
% make sure Y is samples x genes, and is modified by offset

% leave one submission out CV
rng('default');
progfile = 'evaluate_prospective_performance_3_progress_tmp.mat';
if exist(progfile, 'file')
    load(progfile)
    istart = idx + 1;
else
    istart = 1;
    cvi = kfoldcrossvalindbygroup(nfolds, qt.Submission);
    target_proc = zeros(size(Y,1), length(sets));
    target_gene = zeros(size(Y));
    pred_proc = zeros(size(Y,1), length(sets));
    pred_gene = zeros(size(Y));
end

lf = fopen('log.txt', 'w+');
for i = istart : max(cvi)
    fprintf('Fold: %0.0f\n', i);
    fprintf(lf, 'Fold: %0.0f\n', i);
    idx = i; %track the loop variable for progress saving.
    
    % Training
    model = tradict_train_pmvn(Y(cvi ~= i,:), o(cvi ~= i), tids, sets);
    %model = tradict_train( Y(cvi ~= i,:), tids, sets, trainfun, params );
    
    % Build the gene and gene-set target.
    zlag_target = lag_dataset( Y(cvi == i,:), o(cvi == i), 'priors', model.lag_priors); 
    sY_target = standardize(zlag_target, 'mu', model.train_mu, 'std', model.train_sig);
    target_proc(cvi == i, :) = sY_target*model.geneset.coef; % target matrix of pathway/process scores
    target_gene(cvi == i, :) = zlag_target;
    %sY_target = standardize(Y(cvi == i,:), 'mu', model.train_mu, 'std', model.train_sig);
    
    % Formulate the prediction
    [pred_proc(cvi == i ,:), ~, pred_gene(cvi == i, :)] = tradict_predict_pmvn(Y(cvi==i,model.S), o(cvi == i), model);
    
    %pred_proc(cvi == i, :) = predfun(Y(cvi == i,model.S), model.fit);
    
    if mod(i,5) == 0
        save(progfile);
    end
end
fclose(lf);
results.cvi = cvi;
results.target_proc = target_proc;
results.pred_proc = pred_proc;
results.target_gene = target_gene;
results.pred_gene = pred_gene;


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


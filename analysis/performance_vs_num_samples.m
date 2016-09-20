function performance_vs_num_samples( Y, qt, tids, sets, organism )
% Evaluates tradict's performance as a function of increasing number of
% samples.

path(genpath('~/GitHub/tradict'), path)
[ytrain, ytest, ktrain] = partition_data(Y', qt, 0.1);

qt_train = qt(ktrain,:);
qt_test = qt(~ktrain,:);

% Final model for ground truth evaluation
load(['NCBI_SRA_', organism, '_final_tradict_model.mat']);
final_model = model; clear model;


rn = qt_train.release_date_num;
if strcmpi(organism, 'Athaliana'); 
    rn_begin = 734740; %51 samples at this point
else
    rn_begin = 734325;
end
rn_sched = round(logspace(log10(rn_begin), log10(max(rn)), 20)');



% Actual values
T_test = ytest.*repmat(qt_test.spots/1000000,1, size(ytest,2));
o_test = qt_test.spots/1000000;

results = cell(length(rn_sched),1);
for i = 1 : length(rn_sched)
    disp(i);
    smask = rn <= rn_sched(i);
    nsamples(i) = sum(smask);
    nsubs(i) = length(unique( qt_train.Submission(smask)));
    
    ytrain_sub = ytrain(smask,:);
    qtrain_sub = qt_train(smask,:);
    
    T = (ytrain_sub).*repmat(qtrain_sub.spots/1000000,1, size(ytrain_sub,2) );
    o = qtrain_sub.spots/1000000;
    
    % Train Tradict
    model = tradict_train( T, o, tids, sets );
    
    % Formulate prediction
    T_m = T_test(:,model.S);
    pred = tradict_predict( T_m, o_test, model, 'sample_posterior', false );
    
    
    % Actual values
    z = lag_dataset(T_test, o_test, 'priors', final_model.lag_priors);
    zs = standardize(z, 'mu', final_model.train_mu, 'std', final_model.train_sig);
    s = zs*final_model.geneset.coef;
    
    
    
    res.z = z;
    res.s = s;
    res.z_hat = pred.genes.z_hat;
    res.s_hat = pred.programs.s_hat;
    
    
    results{i} = res;
end

save('perf_vs_num_samples_results.mat', 'results', 'ktrain', 'nsamples', 'nsubs', 'rn_sched', '-v7.3');



end


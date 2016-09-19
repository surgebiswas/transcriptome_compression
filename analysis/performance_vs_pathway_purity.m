function performance_vs_pathway_purity( Y, qt, tids, sets )


path(genpath('~/GitHub/tradict'), path)
[ytrain, ytest, ktrain] = partition_data(Y', qt, 0.1);

% Final model for ground truth evaluation
load('NCBI_SRA_Athaliana_final_tradict_model.mat');
final_model = model; clear model;



qt_train = qt(ktrain,:);
qt_test = qt(~ktrain,:);

pswap = fliplr([0 1 2 5 10 20 50 80 100]);


T = (ytrain).*repmat(qt_train.spots/1000000,1, size(ytrain,2) );
o = qt_train.spots/1000000;

T_test = ytest.*repmat(qt_test.spots/1000000,1, size(ytest,2));
o_test = qt_test.spots/1000000;

results = cell(length(pswap),1);
for i = 1 : length(pswap)
    fprintf('Iteration: %0.0f\n', i);
    sets_swapped = perturb_sets(sets,tids,pswap(i));
    
    % Training
    model = tradict_train( T, o, tids, sets_swapped );
    
    % Prediction
    T_m = T_test(:,model.S);
    pred = tradict_predict( T_m, o_test, model, 'calc_credible_intervals', false );
    
    
    % Actual values
    z = lag_dataset(T_test, o_test, 'priors', final_model.lag_priors);
    zs = standardize(z, 'mu', final_model.train_mu, 'std', final_model.train_sig);
    s = zs*final_model.geneset.coef;
    
    res.z = z;
    res.s = s;
    res.z_hat = pred.genes.z_hat;
    res.s_hat = pred.programs.s_hat;
    res.geneset_pexp = model.geneset.pexp;
    
    results{i} = res;
end

save('perf_vs_program_purity.mat', 'results', 'pswap');


    function s = perturb_sets(s, ts, p)
        if p == 0
            ss = s;
            return
        end
        
        for j = 1 : length(s)
            ntoswap = round(length(s{j})*p/100);
            
            
            tdiff = setdiff(ts, s{j});
            idxtoswapout = randsample(length(s{j}), ntoswap);
            idxtoswapin = randsample(length(tdiff), ntoswap);
            
            s{j}(idxtoswapout) = tdiff(idxtoswapin);
        end
        
    end



end


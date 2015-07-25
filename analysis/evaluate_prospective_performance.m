function evaluate_prospective_performance( lY, qt, saveFile )
% lY = [samples x genes] expression matrix. Should be log(TPM + 0.1).
% qt = [samples x fields] query table.
% saveFile = mat file name to which save progress to.


% CONSTANTS/PARAMETERS
params.SAVEMEMORY = false;
params.LAMBDASCHEDULE = logspace(-3, 2, 100);
params.TESTPROP = 0.1; % Proportion of submissions that should be kept in the test set.
params.PROPUNEXPVAR = 0;
params.MAXFEATS = 100;

% 1. Split the data into training and test sets.
[ytrain, ytest, trainind] = partition_data(lY, qt, params.TESTPROP);



% 2. Set up to perform decomposition, tratraining, and tradiction.
% First check to see if the decomposition and tratraining have been
% performed. If so, we can skip them and save some time.
saveFileExists = exist(saveFile, 'file');
if saveFileExists
    load(saveFile);
    if exist('somp', 'var');
        decompositionNeeded = false;
    else
        decompositionNeeded = true;
    end
    
    if exist('model', 'var')
        tratrainingNeeded = false;
    else
        tratrainingNeeded = true;
    end
    
else
    decompositionNeeded = true;
    tratrainingNeeded = true;
end

% 3. Perform training, which includes decomposition and tratraining, as
% necessary.
if decompositionNeeded 
    fprintf('Performing OMP decomposition ... ');
    somp = marker_OMP(standardize(ytrain), params.PROPUNEXPVAR, 'maxfeatures', params.MAXFEATS, ...
        'savememory', params.SAVEMEMORY);
    save(saveFile, 'somp', 'params');
    fprintf('Done.\n');
else
    fprintf('OMP decomposition already performed.\n');
end

if tratrainingNeeded
    fprintf('Performing tratraining ... ');
    model = tratrain(ytrain, ytrain(:,somp.S), 'lambda', params.LAMBDASCHEDULE);
    save(saveFile, 'somp', 'model', 'params');
    fprintf('Done.\n');
else
    fprintf('Tratraining already performed.\n');
end


% 4. Tradict
fprintf('Performing tradiction ... ');
yhat = tradict(ytest(:,somp.S),model);
fprintf('Done.\n');
save(saveFile, 'yhat', 'ytest', 'trainind', 'somp', 'model', 'params');



    

end


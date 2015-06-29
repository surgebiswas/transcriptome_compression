function NCBI_SRA_Athaliana_context_specific_perfomance( lY, qt, latest_training_date )
% How does tradict perform in a context specific setting (e.g. temperature
% perturbed seedlings). In other words while compressing the global
% transcriptome may perform poorly intrasubmission, how do we do if we
% focus e.g. only on seedlings?

% Split into validation and training sets
[lY_validate, lY_train_full, qt_validate] = split_validation_training(lY, qt, latest_training_date);

% Select seedling specific samples.
seedling_mask = select_context_specific_samples(lY_train_full);
lY_train = lY_train_full(seedling_mask,:);


% Marker OMP
PROPVARUNEXPCUTOFF = 0;
MAXFEATS = 100;
mompSaveFile = 'NCBI_SRA_Athaliana_context_specific_performance_MOMP.mat';
if false
    % Marker OMP on seedling specific training data.
    if false
        fprintf('Running Marker OMP on seedling specific training data ...\n');
        somp_seedling = marker_OMP(standardize(lY_train), PROPVARUNEXPCUTOFF, ... 
                    'savememory', true, 'storecrosscorr', true, 'maxfeatures', MAXFEATS);
    end
    
    % Marker OMP on full training data.
    if true
        load(mompSaveFile);
        fprintf('Running Marker OMP on full training data ...\n');
        somp_full = marker_OMP(standardize(lY_train_full), PROPVARUNEXPCUTOFF, ... 
                    'savememory', false, 'storecrosscorr', true, 'maxfeatures', MAXFEATS);
    end
    
    save(mompSaveFile, 'somp_seedling', 'somp_full');
else
    load(mompSaveFile);
end

% Tratrain.
if false
    model_seedling = tratrain(lY_train, lY_train(:,somp_seedling.S), 'lambda', logspace(-3, 2, 20));
    model_full = tratrain(lY_train_full, lY_train_full(:,somp_full.S), 'lambda', logspace(-3, 2, 20));
    save('NCBI_SRA_Athaliana_context_specific_performance_tratrain.mat', 'model_seedling', 'model_full');
else
    load('NCBI_SRA_Athaliana_context_specific_performance_tratrain.mat');
end




% Tradict and Evaluate performance.
seedling_mask_validate = select_context_specific_samples(lY_validate);
lY_validate_seedling = lY_validate(seedling_mask_validate,:);
lY_validate_nonseedling = lY_validate(~seedling_mask_validate,:);
qt_validate_seedling = qt_validate(seedling_mask_validate,:);
qt_validate_nonseedling = qt_validate(~seedling_mask_validate,:);

% Tradict
% Variable names: yhat_<model>_on_<sample_type>
yhat_seedling_on_seedling = tradict(lY_validate_seedling(:,somp_seedling.S), model_seedling);
yhat_seedling_on_nonseedling = tradict(lY_validate_nonseedling(:,somp_seedling.S), model_seedling);
yhat_full_on_seedling = tradict(lY_validate_seedling(:,somp_full.S), model_full);
yhat_full_on_nonseedling = tradict(lY_validate_nonseedling(:,somp_full.S), model_full);


% Performanc evaluation
perf_seedling_on_seedling = evaluate_tradiction(lY_validate_seedling, ...
    yhat_seedling_on_seedling, 'submissionids', qt_validate_seedling.Submission);
perf_seedling_on_nonseedling = evaluate_tradiction(lY_validate_nonseedling, ...
    yhat_seedling_on_nonseedling, 'submissionids', qt_validate_nonseedling.Submission);
perf_full_on_seedling = evaluate_tradiction(lY_validate_seedling, ...
    yhat_full_on_seedling, 'submissionids', qt_validate_seedling.Submission);
perf_full_on_nonseedling = evaluate_tradiction(lY_validate_nonseedling, ...
    yhat_full_on_nonseedling, 'submissionids', qt_validate_nonseedling.Submission);

save('NCBI_SRA_Athaliana_context_specific_performance_statistics.mat', ...
    'perf_seedling_on_seedling', 'perf_seedling_on_nonseedling', ...
    'perf_full_on_seedling', 'perf_full_on_nonseedling');

% Not significant differences between seedling and full model on seedling
% data.

    function [yv, yt, qtv, qtt] = split_validation_training(y, qt, latest_training_date)
        % Create the validation set.
        date_cutoff = datenum(latest_training_date);  
        dn = datenum(qt.release_date);
        q = dn >= date_cutoff;

        yv = y(q,:);
        yt = y(~q,:);
        qtv = qt(q,:);
        qtt = qt(~q,:);
    end

    function mask = select_context_specific_samples(y)
        % Consider updating this to work based on annotation rather than
        % arbitrary PCA cutoffs
        pexp = [];
        coef = [];
        
        load PCA_pexp_vs_eigengene_params;
        s = y*coef(:,1:3);
        
        mask = s(:,1) > -40 & s(:,2) > 20 & s(:,3) < 100;
        
    end




end


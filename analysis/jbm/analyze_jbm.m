clear;
rng('default');

% Initialize environment
completiontext = false;
[~,hn] = unix('hostname'); % Get hostname.
if ~isempty(strfind(hn, '.kd.unc.edu'))
    % We are on killdevil
    homedir = '/proj/dangl_lab/sbiswas';
    completiontext = true;
else 
    % We are working locally
    homedir = '/Users/sbiswas';
    completiontext = false;
end
repo = 'transcriptome_compression/';
datadir = [homedir, '/GitHub/data/', repo, 'jbm/'];
path(genpath([homedir, '/GitHub/', repo]), path);
cd(datadir);


% Set up data.
if false
    load('JBM_assembled_sailfish_output.mat');
    
    % There are some duplicate samples.
    t = tabulate(sids);
    scts = cell2mat(t(:,2));
    torm = steq(sids, t(scts == 2,1));
    
    sids(torm) = [];
    tpm(:,torm) = [];
    
    
    % Collapse isoform table.
    d = mat2dataset(tpm, 'ObsNames', tids, 'VarNames', sids);
    dc = collapse_Athaliana_isoform_table(d);
    
    
    % Keep transcripts used for tradict.
    load('/Users/sbiswas/GitHub/data/transcriptome_compression/Athaliana/NCBI_SRA_Athaliana_full_data_up_to_06Sept2015_quality_filtered.mat');
    dc_sub = dc(tids,:);
    
    
    sids_jbm = get(dc_sub, 'VarNames');
    tids_jbm = get(dc_sub, 'ObsNames');
    lY_jbm = log10(double(dc_sub)' + 0.1);
    
    % Design matrix
    xd = dataset('file', 'JBM_design.txt', 'ReadObsNames', true, 'ReadVarNames', true);
    xd = xd(sids_jbm,:);
    
    % update some of the genotype names. 
    xd.geno(strcmpi(xd.geno, 'oxTCP14-4')) = {'UBQ:TCP14'};
    xd.geno(strcmpi(xd.geno, 'oxTCP14-3')) = {'UBQ:TCP14'};
    xd.geno(strcmpi(xd.geno, 'HopBB1-10')) = {'35S:HopBB1'};
    xd.geno(strcmpi(xd.geno, 'oxJAZ3')) = {'35S:JAZ3'};
    
    
    save('JBM_data_processed.mat', 'sids_jbm', 'tids_jbm', 'lY_jbm', 'xd');
    
    % Some quick exploratory analysis
    s = cmdscale(pdist(standardize(lY_jbm), 'spearman'));
    %[~,s] = pca(standardize(lY_jbm));
    
    fields = {'time', 'treatment', 'rep'};
    colors  = {'r', 'g', 'b', 'k'};
    for i = 1 : length(fields)
        u = unique(eval(['xd.', fields{i}]));
        fprintf('Field: %s\n', fields{i});
        subplot(1,length(fields), i);
        hold on
        for j = 1 : length(u)
            
            if strcmpi(fields{i}, 'rep')
                fprintf('%0.0f -> %s\n', u(j), colors{j});
                m = eval(['xd.', fields{i}]) == u(j);
            else
                fprintf('%s -> %s\n', u{j}, colors{j});
                m = strcmpi(eval(['xd.', fields{i}]), u{j});
            end
            plot(s(m,1), s(m,2), 'o', 'Color', colors{j});
            
            
        end
        title(fields{i});
        if isnumeric(u(1))
            u = cellstr(num2str(u));
        end
        legend(u, 'Location', 'NorthWest');
        fprintf('\n');
    end
    plotSave('figures/exploratory_MDS.png');
    close
    
      
    
    
    
end

load('JBM_data_processed.mat');

% Load/train a tradict model.
if false
    % For now we will use an old set of markers obtained using OMP.
    load('/Users/sbiswas/GitHub/data/transcriptome_compression/Athaliana/NCBI_SRA_Athaliana_full_data_up_to_06Sept2015_quality_filtered.mat');
    load('/Users/sbiswas/GitHub/data/transcriptome_compression/Athaliana/NCBI_SRA_Athaliana_marker_OMP_decomposition_punexp_0.00_maxfeats_500.mat');
    lY = log10(Y' + 0.1);
    
    
    nmarkers = 500;
    markers = somp.S(1:nmarkers);
    
    method = 'ridgefit';
    
    % Loss function
    trimmean_percent = 10;
    lf = @(ytest,yhat) trimmean(abs((ytest(:) - yhat(:))./ytest(:)), ...
        trimmean_percent, 'round', 1); % Robust relative error.
    
    % Cross validation indices.
    % Group by submission.
    nfolds = 10; %length(unique(qt.Submission));
    cvi = kfoldcrossvalindbygroup(nfolds, qt.Submission);
    
    % Representative residual indices.
    rng('default');
    rep_residual = randsample(size(lY,2), round(0.05*size(lY,2)));
    
    if strcmpi(method, 'kernelfit')
        % Kernel fit train parameters
        params_tune{1}{1} = @gaussian_kernel_diagonal; % Kernel
        %params_tune{2} = ([0.1:0.02:0.5]')*sqrt(nmarkers); % loc avg; scale by the dimension of the query space
        params_tune{2} = 2.3*sqrt(nmarkers); % loc reg;
        params_tune{3} = true; % train standardized?
        params_tune{4} = true; % predict original?
        params_tune{5}{1} = 'local_regression'; % ['local_average' | 'local_regression'];

        trainfun = @kernelfit_train;
        predfun = @kernelfit_predict;
    end
    if strcmpi(method, 'ridgefit')
        params_tune{1} = [0, logspace(-6,0,19)]';
        params_tune{2} = true;
        
        % Precompute x'*x and x'*y
        params_tune{3}{1} = ridgefit_precompute;
        params_tune{3}{1}.itr = 1;
        params_tune{3}{1}.precompute = cell(nfolds,2);
        for i = 1 : nfolds
            x = lY(cvi~=i,markers); x = [ones(size(x,1),1),x];
            y = lY(cvi~=i,rep_residual);
            
            params_tune{3}{1}.precompute{i,1} =x'*x;
            params_tune{3}{1}.precompute{i,2} =x'*y;
        end
        
        trainfun = @ridgefit_train;
        predfun = @ridgefit_predict;
    end
    
    
    % Tune the model on a representative residual.
    cvs = cvalidate_tune(lY(:,markers), lY(:,rep_residual), trainfun, ...
        predfun, params_tune, 'crossvalind', cvi, 'loss_fun', lf);
    
    % Build a full model using optimal parameters.
    model = trainfun(lY(:,markers), lY, cvs.param_best_perf);
    
    % Prediction.
    lY_jbm_hat = predfun(lY_jbm(:,markers), model);
    
    save('JBM_tradiction.mat', 'lY_jbm_hat', 'model');
    
end

% Using knowledge based Gene sets.
if true 
    load('JBM_tradiction.mat');
    load('~/Documents/surge/science/gene_ontology/Athaliana_representative_gene_set_01-Mar-2016.mat');
    load('/Users/sbiswas/GitHub/data/transcriptome_compression/Athaliana/NCBI_SRA_Athaliana_gene_set_PC_coefs.mat');
    
    genos = {'Col-0', 'coi1-16', 'npr1-1'}; % ... 
        %'eds16-1', 'jaz3-2', 'UBQ:TCP14', '35S:HopBB1', '35S:JAZ3'};
    
    mask = steq(xd.geno, genos);
    xds = xd(mask,:);

    
    sY_jbm = standardize(lY_jbm(mask,:), 'mu', train_mu, 'std', train_sig);
    sY_jbm_hat = standardize(lY_jbm_hat(mask,:), 'mu', train_mu, 'std', train_sig);
    %[ys, ngenes, engenes] = collapse_to_gene_sets(sY_jbm, tids_jbm, sets);
    %ys_hat = collapse_to_gene_sets(sY_jbm_hat, tids_jbm, sets);
    ys = sY_jbm*gscoef;
    ys_hat = sY_jbm_hat*gscoef;
    
    
    flipv = sign(diag(corr(ys, ys_hat)));
    ys_hat = bsxfun(@times, ys_hat, flipv');
    

    
    o = get(xds, 'ObsNames');
    xs = dataset2cell(xds); xs = xs(2:end, 2:end);
    for i = 1 : size(xs,1)
        if xs{i,end} == 1;
            xs{i,end} = 'A';
        elseif xs{i,end} == 2;
            xs{i,end} = 'B';
        end
    end
    permidx = [2 4 3 1];
    xs = xs(:,permidx);
    vv = {'time', 'rep', 'treat', 'geno'};
    [xss, sidx] = hiersort(xs);
    xssd = cell2dataset(xss, 'ObsNames', o(sidx), 'VarNames', vv(permidx));
    [xbin, vn] = oneHotEncode(xssd, 'intercept', false, 'rmli', false);
    imagesc(1-xbin, [0 1]);colormap(gray);
    plotSave('figures/gene_sets_heatmap_all_major_geno_covars.png');
    close;
    
    
    [ri,ci] = hclust(ys);
    figure
    imagesc(standardize(ys(sidx,ci)), [-3 3]); colormap(prgn);
    axis equal
    axis tight
    plotSave('figures/gene_sets_heatmap_all_major_geno.png');
    close;
    
    figure
    imagesc(standardize(ys_hat(sidx,ci)), [-3 3]); colormap(prgn);
    axis equal
    axis tight
    plotSave('figures/gene_sets_heatmap_col_coi_npr_tradict_estimated.png');
    close;
    
    
end


if false
    load('JBM_tradiction.mat');
    
    t = dataset('file', 'MeJA_up_Col-0_JBM.txt', 'ReadObsNames', false, 'ReadVarNames', false);
    jamarkers = steq(tids_jbm, t.Var1);
    
    t = dataset('file', 'BTH_up_Col-0.txt', 'ReadObsNames', false, 'ReadVarNames', false);
    samarkers = steq(tids_jbm, t.Var1);
    
    
    sja = standardize(lY_jbm(:,jamarkers));
    sja_hat = standardize(lY_jbm_hat(:,jamarkers));
    
    ssa = standardize(lY_jbm(:,samarkers));
    ssa_hat = standardize(lY_jbm_hat(:,samarkers));
    
    % Data point properties
    genos_to_highlight = {'Col-0', 'npr1-1', 'coi1-16', '35S:HopBB1'};
    usecolors = {'k', 'g', 'b', 'r', 'm'};
    mask = steq(xd.geno, genos_to_highlight);
    xdup = xd;
    xdup.geno(~mask) = {'*other'};
    
    
    % JA
    pred_v_actual_density_plot(sja, sja_hat);
    title('MeJA marker prediction accuracy', 'FontSize', 22)
    plotSave('figures/MeJA_marker_prediction_accuracy.png');
    close
    
    figure;
    hold on
    plot([-3 3], [-3 3], '-', 'Color', 0.4*ones(1,3), 'LineWidth', 2)
    [h1,k1] = groupPlot(mean(sja_hat,2), mean(sja,2), xdup, 'usecolors', ...
        usecolors, 'color', 'geno', 'shape', 'treatment', 'markersize', 8);
    
    
    %plot(mean(sja_hat,2), mean(sja,2), 'ok')
    
    text(-1.2, 1.2, sprintf('PCC = %0.2f', corr(mean(sja_hat,2), mean(sja,2))), 'FontSize', 16);
    axis([-1.7 1.7 -1.7 1.7])
    xlabel('Predicted expression', 'FontSize', 18)
    ylabel('Actual expression', 'FontSize', 18)
    set(gca, 'FontSize', 16);
    title('Average(MeJA markers) prediction accuracy', 'FontSize', 20);
    box on;
    axis square;
    plotSave('figures/Average_MeJA_marker_prediction_accuracy.png');
    close
    
    
    % SA
    pred_v_actual_density_plot(ssa, ssa_hat);
    title('BTH marker prediction accuracy', 'FontSize', 22)
    plotSave('figures/BTH_marker_prediction_accuracy.png');
    close
    
    figure;
    hold on
    %plot(mean(ssa_hat,2), mean(ssa,2), 'ok')
    
    plot([-3 3], [-3 3], '-', 'Color', 0.4*ones(1,3), 'LineWidth', 2)
    [h2,k2] = groupPlot(mean(ssa_hat,2), mean(ssa,2), xdup, 'usecolors', ...
        usecolors, 'color', 'geno', 'shape', 'treatment', 'markersize', 8);
    
    text(-1.2, 1.2, sprintf('PCC = %0.2f', corr(mean(ssa_hat,2), mean(ssa,2))), 'FontSize', 16);
    axis([-1.7 1.7 -1.7 1.7])
    xlabel('Predicted expression', 'FontSize', 18)
    ylabel('Actual expression', 'FontSize', 18)
    set(gca, 'FontSize', 16);
    title('Average(BTH markers) prediction accuracy', 'FontSize', 20);
    box on;
    axis square;
    plotSave('figures/Average_BTH_marker_prediction_accuracy.png');
    close
    
    
    
end
% Performs analysis for Parnas et al. 2015. The BDMC CRISPR screeen for Tnf
% induction.

clear;
rng('default');

cd('~/GitHub/data/transcriptome_compression/parnas_bdmc_crispr');

load('SRA248232_assembled_srafish_output.mat');
load('parnas_design.mat');
load('~/GitHub/data/transcriptome_compression/Mmusculus/NCBI_SRA_Mmusculus_full_data_up_to_19Sept2015_quality_filtered.mat', 'tids');



% Pre processing and quality filtering
if false
    % First, get the data and the metadata on the same page.
    Y = s.tpm;
    [all_tids, Y] = collapse_mouse_isoform_table(s.transcript_id, Y);

    qt = read_ncbi_sra_query_table('SRA248232_query_table.csv');
    o = get(xdu, 'ObsNames');
    qts = qt(o,:);
    xdu.depth = qts.spots;
    
    [~, samples_to_keep] = ismember(o, s.ids);
    [~, genes_to_keep] = ismember(tids, all_tids);



    Yt = Y(genes_to_keep, samples_to_keep);
    %ids = s.ids(samples_to_keep);
    xd = xdu;


    % Remove sampels with low average correalation to others
    mr = mean(corr(log10(Yt + 0.1)));
    torm = mr < 0.83;

    Yt(:,torm) = [];
    %sids(torm) = [];
    xd(torm,:) = [];
    Yt = Yt';

    lY = log10(Yt + 0.1);


    % Remove samples with that belong to poorly represented KOs. 
    uk = unique(xd.KO);
    genos_torm = [];
    for i = 1 : length(uk)
        mask = strcmpi(xd.KO, uk{i});

        if length(unique(xd.treatment(mask))) < 4
            genos_torm = [genos_torm; uk(i)];
        end

    end

    torm = steq(xd.KO, genos_torm);
    xd(torm,:) = [];
    lY(torm,:) = [];
    Yt(torm,:) = [];
    
    
    % Add pathway and regulator class information
    canonical = {'Lbp', 'Ly96', 'Tlr4', 'Cd14', 'Tirap', 'Myd88', 'Ticam2', ...
        'Ticam1', 'Irak4', 'Irak1', 'Ripk1', 'Traf6', 'Tab1', 'Tab2', 'Map3k7', ...
        'Irf5', 'Ikbkg', 'Chuk', 'Nfkb1', 'Ikbkb', 'Map3k8', 'Nfkbia', 'Rela', 'Tnf'};
    ost = {'Srp54c', 'Srpr', 'Sec61', 'Rpn2', 'Rpn1', 'Ddost', 'Dad1', 'Alg2', 'Hsp90b1', ...
        'Sec13'};
    paf = {'Paf1', 'Ctr9', 'Wdr61', 'Rtf1', 'Leo1'};
    kc = dataset('file', 'KO_categories_nonredundant.csv', 'ReadObsNames', true, 'ReadVarNames', false, 'Delimiter', ',');
    kc( sum(double(kc),2) == 0, : ) = [];
    ko2direc = containers.Map;
    o = get(kc, 'ObsNames');
    for i = 1 : length(kc)
        if kc.Var2(i) == 1
            ko2direc(o{i}) = 1; % positive regulator
        elseif kc.Var3(i) == 1
            ko2direc(o{i}) = -1; % negative regulator
        end
    end
    
    pathway = cell(size(xd,1),1);
    regclass = cell(size(xd,1),1);
    for i = 1 : length(pathway)
        if any(strcmpi(canonical, xd.KO{i}))
            pathway{i} = 'canonical';
        elseif any(strcmpi(ost, xd.KO{i}))
            pathway{i} = 'ost';
        elseif any(strcmpi(paf, xd.KO{i}))
            pathway{i} = 'paf';
        else
            pathway{i} = 'other';
        end
        
        if ko2direc.isKey(xd.KO{i})
            regclass{i} = ko2direc(xd.KO{i});
        else
            regclass{i} = 0;
        end
    end
    
    toadd = cell2dataset([regclass, pathway], 'VarNames', {'reg_class', 'pathway'}, 'ObsNames', get(xd, 'ObsNames'));
    xd = [xd, toadd];
    
    
    % focus on KOs where we know whether its a + or - regulator and in
    % which Brefeldin A has not been added (since those are the RNA-seq
    % results they describe in the paper). 
    % Also remove negative regulators from the analysis since a saturating
    % dose of LPS given. Also it does not looks like that the negative
    % regulators strongly affect response to LPS or cell response to Tnf.
    torm = (xd.reg_class == 0 & ~strcmpi(xd.KO, 'none')) | ~strcmpi(xd.condition, 'NA') | xd.reg_class == -1;
    xd(torm,:) = [];
    lY(torm,:) = [];
    Yt(torm,:) = [];
    
    xd.treatment = regexprep(xd.treatment, 'LPS  ', 'LPS ');
    
    xd.source_name = regexprep(xd.source_name, 'T\d', '');
    us = unique(xd.source_name);
    
    save('parnas_processed.mat', 'lY', 'Yt', 'xd', 'tids');
else
    load('parnas_processed.mat');
    
end

load('~/GitHub/transcriptome_compression/analysis/gene_ontology/Mmusculus_representative_gene_set_02-Apr-2016.mat');
load('~/GitHub/data/transcriptome_compression/Mmusculus/NCBI_SRA_Mmusculus_final_tradict_model_marker_genes_expression_mat.mat');
load('~/GitHub/data/transcriptome_compression/Mmusculus/NCBI_SRA_Mmmusculus_final_tradict_model.mat');


splits = regexpi(tids, ',', 'split');
sp = reshape([splits{:}], 2, length(tids))';
stids = tids;
tids = sp(:,2);

% Create the design matrix for batch effect correction.
if true % leave as true
    % Estimate major effects
    [~, ss] = pca(lY, 'NumComponents', 10);
    batch1 = ss(:,1) < -20;
    batch2 = ss(:,1) > -20 & ss(:,3) < -4;
    batch3 = ss(:,1) > -20 & ss(:,3) > -4 & ss(:,3) < 7.2;
    batch4 = ss(:,1) > -20 & ss(:,3) > 7.2;
    treat = zeros(size(lY,1),1);
    treat(strcmpi(xd.treatment, 'LPS 2 hours')) = 1;
    treat(strcmpi(xd.treatment, 'LPS 4 hours')) = 2;
    treat(strcmpi(xd.treatment, 'LPS 6 hours')) = 3;
    
    %xeff = [ones(length(batch1),1), batch1, batch2, batch3, treat];
    xeff = [ones(length(batch1),1), batch1, treat];
end


% Form the predictions and lag the actual.
load('~/GitHub/data/transcriptome_compression/Mmusculus/NCBI_SRA_Mmusculus_final_tradict_model.mat');
if false

    t_m = Yt(:, model.S).*repmat(xd.depth/1000000,1, length(model.S));
    o = xd.depth/1000000;

    [ s_hat, Yt_hat, z_hat ] = tradict_predict_pmvn( t_m, o, model );
    [ s_hat_d, Yt_hat_d, z_hat_d] = tradict_predict_pmvn( t_m, o, model, 'learn_latent_diag', true);
    save('parnas_tradiction.mat', 's_hat', 'Yt_hat', 'z_hat');
end

% Do comparison 
load('parnas_tradiction.mat')

% Create the actual gene-set scores
if true
    t = Yt.*repmat(xd.depth/1000000,1, size(Yt,2));
    o = xd.depth/1000000;
    
    z = lag_dataset(t, o, 'priors', model.lag_priors);
    zs = standardize(z, 'mu', model.train_mu, 'std', model.train_sig);
    s = zs*model.geneset.coef;
    
    
    
end






% Gene set adjustments 
sa = parnas_batch_correct(xeff, s, 2);
s_hata = parnas_batch_correct(xeff, s_hat, 2);
%s_hata_d = parnas_batch_correct(xeff, s_hat_d, 2);

sas = standardize(sa);
s_hatas = standardize(s_hata);
%s_hatas_d = standardize(s_hata_d);

% Gene adjustments
za = parnas_batch_correct(xeff, z, 2);
z_hata = parnas_batch_correct(xeff, z_hat, 2);



% plot for response to LPS
if true
    cm = lines;
    ut = unique(xd.treatment); ut = [ut(4); ut(1:3)];
    lpsi = find(strcmpi(setnames(:,1), 'cellular response to lipopolysaccharide'));
    
    hold on
    for i = 1 :length(ut)
        m = strcmpi(xd.treatment, ut{i});
        l(i) = plot(s_hatas(m,lpsi), sas(m,lpsi), 'o', 'MarkerFaceColor', cm(i,:), 'Color', cm(i,:));
        disp(corr(sas(m,lpsi), s_hatas(m,lpsi), 'type', 'spearman'));
    end
    
    sf = get_standard_figure_font_sizes;
    axis square;
    legend(l, ut, 'Location', 'SouthEast');
    axis([-2.5 2.5 -2.5 2.5]);
    set(gca, 'FontSize', sf.axis_tick_labels);
    xlabel('Predicted score', 'FontSize', sf.axis_labels);
    ylabel('Actual score', 'FontSize', sf.axis_labels);
    plotSave('LPS_pathway_score.png'); close;
    
    
    
    
    
end



% Perform TPM correction for low read counts
% t = 10.^lY(:,model.S) - 0.1;
% 
% mu = mean(lY(:,model.S));
% Sigma = cov(lY(:,model.S)); %cov(log(ysub + 0.1));
% zhat = learn_pmvn_z(t, mu, Sigma, ones(size(t,1),1));
% 
% % Sigma = cov(zhat);
% % mu = mean(zhat);
% % % zhat = learn_pmvn_z(t, mu, Sigma, ones(size(t,1),1));
% % 
% % lY_c = log10(exp(zhat) + 0.1);
% 
% % Gene set adjustments - predicted
% gYh = ridgefit_predict(lY(:,model.S), model.fit); % ow lY(:,model.S)
% b = estimate_effects(xeff, gYh);
% gYhadj = gYh - xeff(:,2:4)*b(2:4,:);
% syha = standardize(gYhadj);
% 
% % Gene set adjustments - actual
% lYh = ridgefit_predict(lY(:,model.S), model.fit_allgenes); % ow lY(:,model.S)
% b = estimate_effects(xeff, lYh);
% lYhadj = lYh - xeff(:,2:4)*b(2:4,:);
% lyha = standardize(lYhadj);



if false
    ir = {
    'immune response';
    'neutrophil chemotaxis';
    'chemotaxis';
    'positive regulation of ERK1 and ERK2 cascade';
    'defense response';
    'cellular response to tumor necrosis factor';
    'positive regulation of inflammatory response';
    'positive regulation of T cell proliferation';
    'regulation of inflammatory response';
    'inflammatory response';
    'cellular response to interferon-gamma';
    'negative regulation of inflammatory response';
    'defense response to Gram-positive bacterium';
    'response to lipopolysaccharide';
    'defense response to bacterium';
    'positive regulation of interferon-gamma production';
    'positive regulation of JNK cascade';
    'positive regulation of tumor necrosis factor production';
    'immune system process';
    'activation of MAPK activity';
    'cytokine-mediated signaling pathway';
    'cellular response to lipopolysaccharide';
    'innate immune response';
    'positive regulation of interleukin-6 production';
    'positive regulation of NF-kappaB transcription factor activity',
    'response to virus'};

    [~,proc2keep] = ismember(ir,setnames(:,1));

    m = steq(setnames(:,1), ir);
    imgenes = sum(model.geneset.coef(:,m),2) ~= 0;
    [kop,ci] = plot_immune_heatmaps(xd, standardize(lY(:,imgenes)), snames, [], []);
    
    
    
    
    
    
end

% Differential pathway expression analysis.
if false
    % Build design matrix.
    time = zeros(size(xd,1),1);
    time(strcmpi(xd.treatment, 'No LPS')) = 0; 
    time(strcmpi(xd.treatment, 'LPS 2 hours')) = 2; 
    time(strcmpi(xd.treatment, 'LPS 4 hours')) = 4; 
    time(strcmpi(xd.treatment, 'LPS 6 hours')) = 6; 
    
    rc = xd.reg_class;
    
    X = [time, rc];
    
    % Build the linear models and obtain p-values]
    % This is a well behaved linear model with two terms and low
    % multicolinearity
    pval_true = zeros(size(X,2) + 1, size(sas,2));
    pval_pred = pval_true;
    for i = 1 : size(sas,2)
        stats_true = regstats(sas(:,i), X, 'linear');
        stats_pred = regstats(s_hatas(:,i), X, 'linear');
        
        pval_true(:,i) = stats_true.tstat.pval;
        pval_pred(:,i) = stats_pred.tstat.pval;
        disp(i);
    end
    
    % Adjust to FDR for reg class coefficint 
    fdr_true = mafdr(pval_true(3,:)', 'BHFDR', true);
    fdr_pred = mafdr(pval_pred(3,:)', 'BHFDR', true);
    
    % Make a report.
    [~,sidx] = sort(fdr_pred);
    
    mr = [num2cell([fdr_true(sidx), fdr_pred(sidx)]), setnames(sidx,1)];
    
    
    labels = fdr_true < 0.01;
    scores = 1 - fdr_pred;
    [xx,yy, th, auc] = perfcurve(labels, scores, 1);
    [~,mind] = min(abs(th - (1 - 0.01)))
    
    % Figures
    figure;
    plot(xx, yy, '-k', 'LineWidth', 3);
    hold on
    plot([xx(mind), xx(mind)], [0 yy(mind)], '--r', 'LineWidth', 3);
    plot([0, xx(mind)], [yy(mind) yy(mind)], '--r', 'LineWidth', 3);
    plot([xx(mind)], [yy(mind)], 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 12);
    axis square;

    sf = get_standard_figure_font_sizes;
    axis([0 1 0 1]);
    set(gca, 'FontSize', sf.axis_tick_labels);
    xlabel('Estimated FPR', 'FontSize', sf.axis_labels);
    ylabel('Estimated TPR', 'FontSize', sf.axis_labels);
    plotSave('roc_curve.png');
    close
    
    cm = parula;
    v = venn([sum(fdr_true < 0.01), sum(fdr_pred < 0.01)], sum(fdr_pred < 0.01 & fdr_true < 0.01));
    set(v(1), 'FaceColor', cm(10,:))
    set(v(2), 'FaceColor', cm(end-10,:))
    axis equal
    axis off
    plotSave('venn.png');
    close
    
    
    % heatmaps
    mask = fdr_pred < 0.01 | fdr_true < 0.01;
    [pris, ci] = plot_immune_heatmaps(xd, sas(:,mask));
    plotSave('DE_genesets_heatmap_true.png');
    close
    
    plot_immune_heatmaps(xd, s_hatas(:,mask), 'pris', pris, 'ci', ci);
    plotSave('DE_genesets_heatmap_predicted.png');
    close
    
    dek = [fdr_pred < 0.01 , fdr_true < 0.01];
    dek = dek(mask,:);
    imagesc(1 - dek(ci,:), [0 1]); colormap(gray);
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca, 'TickLength', [0 0]);
    daspect([1 4 1]);
    plotSave('DE_annotations.png');
    close;
end







% scratch
if false
    [~,s,~,~,pexp] = pca(lY, 'NumComponents', 10);


    v = xd.treatment;
    ut = unique(v);
    figure;
    hold on
    for i = 1 : length(ut)
        mask = strcmpi(v, ut{i}); 
        plot3(s(mask,1), s(mask,2), s(mask,3), '.', 'Color', rand(1,3));
    end

end
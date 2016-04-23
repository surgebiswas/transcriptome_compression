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

    o = get(xdu, 'ObsNames');
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

    lY = log10(Yt' + 0.1);


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
    
    
    save('parnas_processed.mat', 'lY', 'xd', 'tids');
else
    load('parnas_processed.mat');
    
end

load('~/GitHub/transcriptome_compression/analysis/gene_ontology/Mmusculus_representative_gene_set_02-Apr-2016.mat');
load('~/GitHub/data/transcriptome_compression/Mmusculus/NCBI_SRA_Mmmusculus_final_tradict_model.mat');


% Estimate major effects
[~, s] = pca(lY, 'NumComponents', 10);
batch1 = s(:,1) < -20;
batch2 = s(:,1) > -20 & s(:,3) < -4;
batch3 = s(:,1) > -20 & s(:,3) > -4 & s(:,3) < 7.2;
batch4 = s(:,1) > -20 & s(:,3) > 7.2;
treat = zeros(size(lY,1),1);
treat(strcmpi(xd.treatment, 'LPS  2 hours')) = 1;
treat(strcmpi(xd.treatment, 'LPS  4 hours')) = 2;
treat(strcmpi(xd.treatment, 'LPS  6 hours')) = 3;

xeff = [ones(length(batch1),1), batch1, batch2, batch3, treat];



gY = lY*model.geneset.coef;
b = estimate_effects(xeff, gY);
gYadj = gY - xeff(:,2:4)*b(2:4,:);
sya = standardize(gYadj);

gYh = ridgefit_predict(lY(:,model.S), model.fit);
b = estimate_effects(xeff, gYh);
gYhadj = gYh - xeff(:,2:4)*b(2:4,:);
syha = standardize(gYhadj);

% Differential pathway expression analysis.
if true
    % Build design matrix.
    time = zeros(size(xd,1),1);
    time(strcmpi(xd.treatment, 'No LPS')) = 0; 
    time(strcmpi(xd.treatment, 'LPS  2 hours')) = 2; 
    time(strcmpi(xd.treatment, 'LPS  4 hours')) = 4; 
    time(strcmpi(xd.treatment, 'LPS  6 hours')) = 6; 
    
    rc = xd.reg_class;
    
    X = [time, rc];
    
    % Build the linear models and obtain p-values]
    % This is a well behaved linear model with two terms and low
    % multicolinearity
    pval_true = zeros(size(X,2) + 1, size(gYadj,2));
    pval_pred = pval_true;
    for i = 1 : size(gYadj,2)
        stats_true = regstats(gYadj(:,i), X, 'linear');
        stats_pred = regstats(gYhadj(:,i), X, 'linear');
        
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
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
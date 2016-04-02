clear;
rng('default');
cd('/Users/sbiswas/Documents/surge/science/gene_ontology');

% Defining a representative set of BIOLOGICAL PROCESSES GO categories as
% gene sets for A. thaliana. 


load Athaliana_gene_association_16Feb2016.mat

GOdbpath = '/Users/sbiswas/Documents/surge/science/gene_ontology/gene_ontology_db_17-Feb-2016.obo';
GO = geneont('File', GOdbpath);

genes = unique({ann.DB_Object_Name}');
goids = unique([ann.GOid]);

% fast indexing
genes2idx = containers.Map;
idx2genes = containers.Map('KeyType', 'double', 'ValueType', 'char');
for i = 1 : length(genes)
    genes2idx(genes{i}) = i; 
    idx2genes(i) = genes{i};
end

go2idx = containers.Map('KeyType', 'double', 'ValueType', 'double');
idx2go = containers.Map('KeyType', 'double', 'ValueType', 'double');
for i = 1 : length(goids)
    go2idx(goids(i)) = i; 
    idx2go(i) = goids(i);
end


ggt = zeros( numel(genes), numel(goids) );
for i = 1 : numel(ann)
    rowidx = genes2idx(ann(i).DB_Object_Name);
    colidx = go2idx(ann(i).GOid);
    ggt(rowidx,colidx) = ggt(rowidx,colidx) + 1;
end
sg = sum(ggt > 0);

if true
    % Range analysis. Want to maximize number of genes covered and minimize of go
    % terms. 
    minrange = [1 10 50 100];
    maxrange = [500 1000 1500 2000 5000 10000];
    for i = 1 : length(minrange)
        for j = 1 : length(maxrange)
            gr = sg >= minrange(i) & sg <= maxrange(j);

            genescovered = nnz(sum(ggt(:,gr),2));
            numgoterms = nnz(gr);

            fprintf('Min GO Term size: %0.0f\tMax GO Term Size: %0.0f\tNum. Genes Covered: %0.0f\tNum. GO Terms: %0.0f\n', ...
                minrange(i), maxrange(j), genescovered, numgoterms);
        end
    end
    % seems min size = 50, max size = 2000 is the best. 
    % Min GO Term size: 50	Max GO Term Size: 2000	Num. Genes Covered: 15269	Num. GO Terms: 168
end

% Representative GO term set.
minfinal = 50;
maxfinal = 2000;
grfinal = sg >= minfinal & sg <= maxfinal;

goids_final = [];
grfinalidx = find(grfinal);
for i = 1 : length(grfinalidx)
    goids_final = [goids_final, idx2go(grfinalidx(i))];
end

ggt_final = ggt(:,grfinal);
ggt_binary = ggt_final > 0;

setnames = cell(length(goids_final),2);
sets = cell(length(goids_final),1);
for i = 1 : length(goids_final)
    setnames(i,:) = {GO(goids_final(i)).Terms.name, GO(goids_final(i)).Terms.definition};
    sets{i} = genes(ggt_binary(:,i));
end

save(['Athaliana_representative_gene_set_', datestr(today), '.mat'], 'sets', 'setnames', 'goids_final');





% uids = unique([ann.GOid]);
% names = cell(numel(uids),1);
% for i = 1 : length(names)
%     names{i} = GO(uids(i)).Terms.name;
%     if mod(i,500) == 0 
%         disp(i)
%     end
% end
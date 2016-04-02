function build_association_mat_file( aspect, organism )

switch organism
    case 'Arabidopsis thaliana'
        
        association_file = '/Users/sbiswas/Documents/surge/science/gene_ontology/gene_association_16Feb2016.tair';
        gene_id_field = 'DB_Object_Name';
        matfile = '/Users/sbiswas/Documents/surge/science/gene_ontology/Athaliana_gene_association_16Feb2016.mat';
        
        fprintf('Loading association file ... ');
        ann = goannotread(association_file,'Aspect', aspect,...
                             'Fields', {gene_id_field,'GOid'});
        fprintf('Done.\n');

        % amap is a hash that maps gene IDs (with at least 1 GO annotation) to 
        % their GO IDs. This is built using the assocation file.
        fprintf('Building gene ID -> GO ID hash ... ');
        amap = containers.Map();
        for i = 1:numel(ann)
            key = eval(['ann(i).', gene_id_field]);
            if isKey(amap,key)
                amap(key) = [amap(key) ann(i).GOid];
            else
                amap(key) = ann(i).GOid;
            end
        end
        fprintf('Done.\n');

        save(matfile, 'ann', 'amap');

        
        
        
    case 'Mus musculus'
        
        association_file = '/Users/sbiswas/Documents/surge/science/gene_ontology/gene_association_29Mar2016.mgi';
        gene_id_field = 'DB_Object_Symbol';
        matfile = '/Users/sbiswas/Documents/surge/science/gene_ontology/Mmusculus_gene_association_29Mar2016.mat';
        
        fprintf('Loading association file ... ');
        ann = goannotread(association_file,'Aspect', aspect,...
                             'Fields', {gene_id_field,'GOid'});
        fprintf('Done.\n');

        % amap is a hash that maps gene IDs (with at least 1 GO annotation) to 
        % their GO IDs. This is built using the assocation file.
        fprintf('Building gene ID -> GO ID hash ... ');
        amap = containers.Map();
        for i = 1:numel(ann)
            key = eval(['ann(i).', gene_id_field]);
            if isKey(amap,key)
                amap(key) = [amap(key) ann(i).GOid];
            else
                amap(key) = ann(i).GOid;
            end
        end
        fprintf('Done.\n');

        save(matfile, 'ann', 'amap');
       

end


end


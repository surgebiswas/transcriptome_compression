function s = assemble_srafish_output( srafish_outdir, sra_ids, varargin )
% srafish_outdir = output directory of srafish
% sra_ids = the SRA run ID's which should be included in the assembly.

quickread = setParam(varargin, 'quickread', true);

old = cd(srafish_outdir);

s.ids = sra_ids; if ~isrow(s.ids); s.ids = s.ids'; end
s.mapped_ratio = zeros(1,length(sra_ids));
tic;
BATCHSIZE = 1;
for i = 1 : length(sra_ids)
    cd(sra_ids{i});

    

    if i == 1
        % Read in the full table the first pass so that we can get the
        % transcript IDs. 
        q = read_sailfish_output('quant_bias_corrected.sf', 'quickread', false);
        
        s.tpm = [q.table.TPM, zeros(length(q.table.TPM),length(sra_ids)-1)];
        s.estNumReads = [q.table.EstimatedNumReads, zeros(length(q.table.TPM),length(sra_ids)-1)];
        s.transcript_id = get(q.table, 'ObsNames');
    else
        % In all subsequnt passes just read the TPM column.
        % Read quickly if user requests.
        q = read_sailfish_output('quant_bias_corrected.sf', 'quickread', quickread);
        
        s.tpm(:,i) = q.table(:,1);
        s.estNumReads(:,i) = q.table(:,2);
    end
  
    rci = dataset('file', 'reads.count_info', 'ReadObsNames', true, 'ReadVarNames', false);
    s.mapped_ratio(i) = rci.Var2(4);
  
    cd('..');
  
  
    if mod(i,BATCHSIZE) == 0
        t = toc;
        fprintf('%0.0f record(s) assembled. Last %0.0f records took %0.2f seconds to assemble.\n', i, BATCHSIZE, t);
        tic;
    end
end

s.depth = sum(s.estNumReads);



cd(old)



end
function s = assemble_srafish_output( srafish_outdir, sra_ids )
% srafish_outdir = output directory of srafish
% sra_ids = the SRA run ID's which should be included in the assembly.

old = cd(srafish_outdir);

s.ids = sra_ids; if ~isrow(s.ids); s.ids = s.ids'; end
s.mapped_ratio = zeros(1,length(sra_ids));
tic;
BATCHSIZE = 1;
for i = 1 : length(sra_ids)
    cd(sra_ids{i});

    q = read_sailfish_output('quant_bias_corrected.sf');

    if i == 1
        s.tpm = [q.table.TPM, zeros(length(q.table.TPM),length(sra_ids)-1)];
        s.rpkm = [q.table.RPKM, zeros(length(q.table.RPKM),length(sra_ids)-1)];
        s.estNumReads = [q.table.EstimatedNumReads, zeros(length(q.table.EstimatedNumReads),length(sra_ids)-1)];
        s.transcript_id = get(q.table, 'ObsNames');
    else
        s.tpm(:,i) = q.table.TPM;
        s.rpkm(:,i) = q.table.RPKM;
        s.estNumReads(:,i) = q.table.EstimatedNumReads;
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
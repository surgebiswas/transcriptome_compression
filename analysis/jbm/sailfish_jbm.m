
old = cd('/proj/dangl_lab/sbiswas/GitHub/data/transcriptome_compression/jbm/count_files');
d = dataset('file', 'jbm_to_sailfish.txt', 'ReadVarNames', false, 'ReadObsNames', false);

for i = 1 : length(d.Var1{i})
    % Copy the file
    cmd = sprintf('cp %s .', d.Var1{i});
    fprintf('%s\n', cmd);
    %system(cmd);
    
    [p,f,e] = fileparts(d.Var1{i});
    fname = strrep(d.Var1{i}, [p, '/'], '');
    fnamefasta = strrep(fname, '.gz', '');
    
    % Gunzip
    cmd = sprintf('gunzip %s', fname);
    fprintf('%s\n', cmd);
    %system(cmd);
    
    
    % Sailfish
    splits = regxpi(fname, '_', 'split');
    outdir = [splits{1}, '_', splits{2}];
    mkdir(outdir);
    cmd = sprintf('sailfish quant -i %s -l ''T=SE:S=U'' -r %s -o %s -p 6 &> %s', ...
        '/proj/dangl_lab/sbiswas/sailfish_indexes/tair10', fnamefasta, outdir, [outdir, '/sf.out']);
    bcmd = sprintf('bsub -q day -M 30 -o %s -e %s -n 6 -R "span[hosts=1]" %s', ...
        [outdir, '/lsf.out'], [outdir, '/lsf.err'], cmd);
    fprintf('%s\n', cmd);
    %system(cmd);
end

cd(old);
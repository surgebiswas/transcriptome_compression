
old = cd('/proj/dangl_lab/sbiswas/GitHub/data/transcriptome_compression/jbm/count_files');
d = dataset('file', 'jbm_to_sailfish.txt', 'ReadVarNames', false, 'ReadObsNames', false);
cmdfile = 'jbm_sailfish_cmd_calls.txt';
fid = fopen(cmdfile, 'w+');

for i = 1 : length(d.Var1)
    % Copy the file
    cmd = sprintf('cp %s .', d.Var1{i});
    fprintf('%s\n', cmd);
    system(cmd);
    
    [p,f,e] = fileparts(d.Var1{i});
    fname = strrep(d.Var1{i}, [p, '/'], '');
    fnamefasta = strrep(fname, '.gz', '');
    
    % Gunzip
    cmd = sprintf('gunzip %s', fname);
    fprintf('%s\n', cmd);
    system(cmd);
    
    % Make and move file to output directory.
    splits = regexpi(fname, '_', 'split');
    outdir = [splits{1}, '_', splits{2}];
    mkdir(outdir);
    cmd = sprintf('mv %s %s', fnamefasta, outdir);
    fprintf('%s\n', cmd);
    system(cmd);
    
    
    % Sailfish
    cmd = sprintf('sailfish quant -i %s -l ''T=SE:S=U'' -r %s -o %s -p 6 &> %s', ...
        '/proj/dangl_lab/sbiswas/sailfish_indexes/tair10', [outdir, '/', fnamefasta], outdir, [outdir, '/sf.out']);
    bcmd = sprintf('bsub -q day -M 30 -o %s -e %s -n 6 -R "span[hosts=1]" "%s"', ...
        [outdir, '/lsf.out'], [outdir, '/lsf.err'], cmd);
    fprintf('%s\n', bcmd);
    %system(cmd);
    
    fprintf(fid, '%s\n', bcmd);

end
fclose(fid);
cd(old);
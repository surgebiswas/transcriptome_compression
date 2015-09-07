function cmd = create_srafish_calls_on_subtables(qt_subtable_direc, organism, varargin)

    cd(qt_subtable_direc);
    %Assumes every file in this directory is query tables.
    queue = setParam(varargin, 'queue', 'day');
    mem = setParam(varargin, 'memory', 30);
    nt = setParam(varargin, 'numthreads', 6);
    donotprocess = setParam(varargin, 'donotprocess', []);
    opssh = setParam(varargin, 'openssh', ...
        '/proj/dangl_lab/sbiswas/GitHub/transcriptome_compression/data_download/asperaweb_id_dsa.openssh');
    minreads = setParam(varargin, 'minreads', 4000000);
    maxreads = setParam(varargin, 'maxreads', 70000000);
    
    if strcmpi(organism, 'Mmusculus')
        index = '/proj/dangl_lab/sbiswas/sailfish_indexes/gencode_vM5';
        outdir = '/lustre/scr/s/b/sbiswas/NCBI_SRA/Mmusculus/';
    elseif strcmpi(organism, 'Athaliana')
        index = '/proj/dangl_lab/sbiswas/sailfish_indexes/tair10';
        outdir = '/lustre/scr/s/b/sbiswas/NCBI_SRA/Athaliana/';
    end

    td = [pwd, '/'];
    listing = dir(pwd);
    listing(1:2) = [];
    cmd = cell(length(listing),1);
    for i = 1 : length(listing)
        f = listing(i).name;

        bsubcall = sprintf('bsub -q %s -M %0.0f -o %s -e %s -n %0.0f -R "span[hosts=1]"', ...
                    queue, mem, [outdir, f, '.out'], [outdir, f, '.err'],  nt);
        
        if isempty(donotprocess)
            cmd{i} = sprintf('%s "~/GitHub/transcriptome_compression/data_download/srafish.pl -t %s -i %s -s %s -n %0.0f -o %s -m %0.0f -a %0.0f"', ...
                bsubcall, [td,f], index, opssh, nt, outdir, minreads, maxreads);
        else
            cmd{i} = sprintf('%s "~/GitHub/transcriptome_compression/data_download/srafish.pl -t %s -i %s -s %s -n %0.0f -o %s -m %0.0f -a %0.0f -d %s"', ...
                bsubcall, [td,f], index, opssh, nt, outdir, minreads, maxreads, donotprocess);
        end

    end
    
end

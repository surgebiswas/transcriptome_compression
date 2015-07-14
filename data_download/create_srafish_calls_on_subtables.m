function cmd = create_srafish_calls_on_subtables(qt_subtable_direc)

    cd(qt_subtable_direc);
    %Assumes every file in this directory is query tables.
    queue = 'day';
    mem = 30;
    nt = 6;
    index = '/proj/dangl_lab/sbiswas/transcriptome_compression/NCBI_SRA/sailfish_indexes/gencode_vM5';
    opssh = '/nas02/home/s/b/sbiswas/GitHub/transcriptome_compression/data_download/asperaweb_id_dsa.openssh';
    donotprocess = '/proj/dangl_lab/sbiswas/transcriptome_compression/NCBI_SRA/Mmusculus/download_04June2015/successfully_processed_list.txt';
    outdir = '/lustre/scr/s/b/sbiswas/NCBI_SRA/Mmusculus/';
    minreads = 4000000;
    maxreads = 70000000;


    td = [pwd, '/'];
    listing = dir(pwd);
    listing(1:2) = [];
    cmd = cell(length(listing),1);
    for i = 1 : length(listing)
        f = listing(i).name;

        bsubcall = sprintf('bsub -q %s -M %0.0f -o %s -e %s -n %0.0f -R "span[hosts=1]"', ...
                    queue, mem, [outdir, f, '.out'], [outdir, f, '.err'],  nt);

        cmd{i} = sprintf('%s "~/GitHub/transcriptome_compression/data_download/srafish.pl -t %s -i %s -s %s -n %0.0f -o %s -m %0.0f -a %0.0f -d %s"', ...
                bsubcall, [td,f], index, opssh, nt, outdir, minreads, maxreads, donotprocess);

    end
    
end

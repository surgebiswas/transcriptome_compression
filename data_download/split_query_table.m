function split_query_table( qtfile, outdir, nsplit )
% Splits the query table specified by qtfile into multiple subtables. 
% Subtables are N.csv where N is the
% subtable index.
% nsplit = number of queries to include per sub table.
old = cd(outdir);

qt = read_ncbi_sra_query_table(qtfile);
fid = fopen(qtfile);
hline = fgetl(fid); % header line

lines = cell(size(qt,1),1);
for i = 1 : length(lines)
    lines{i} = fgetl(fid);
end
fclose(fid);

% Randomize the lines so that if patches of bad samples get mixed up.
rlines = lines(randperm(length(lines)));

% Print the subtables. 
nfout = 0;
tname = tempname;
fout = fopen(tname, 'w+');
for i = 1 : length(rlines)
    q = ceil(i/nsplit);
    if q ~= nfout
        nfout = q;
        fclose(fout); % close previous handle
        fout = fopen([num2str(nfout), '.csv'], 'w+'); % open new handle
        fprintf(fout, '%s\n', hline); % print header line 
    end
    
    fprintf(fout, '%s\n', rlines{i});
end
system(['rm ', tname]);
fclose(fout);

cd(old);

end


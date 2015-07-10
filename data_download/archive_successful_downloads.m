function archive_successful_downloads( processed_file, archive_basename )
% Archives successful downloads by srafish.
%
% processed_file - contains SRA IDs of processed entries on separate lines.
% archive_basename - file basename of .tgz archive. Final archive will be
%                    [basename.tgz]

fid = fopen(processed_file);
t = textscan(fid, '%s\n');
fclose(fid);

ids = t{1};

mkdir(archive_basename)
fprintf('Moving data ... ');
for i = 1 : length(ids);
    system(sprintf('mv %s %s', ids{i}, archive_basename));
end
fprintf('Done.\n');

fprintf('Compressing ... ');
system(sprintf('tar -zcvf %s.tgz %s', archive_basename, archive_basename)); 
fprintf('Done.\n');




end


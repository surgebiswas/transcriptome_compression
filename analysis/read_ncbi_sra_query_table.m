function [ qt ] = read_ncbi_sra_query_table( qtfile )

[~, result] = system( ['wc -l ', qtfile] );
numlines = str2double( strrep(result, qtfile, '') );



fid = fopen(qtfile);
l = fgetl(fid);

headers = regexp(l, ',', 'split');
NCOLRM = 3;
headers(1:NCOLRM) = []; % remove run id, release, and load date.

l = fgetl(fid);
release_dates = cell(numlines-1,1);
load_dates = cell(numlines-1,1);
md = cell(numlines-1,length(headers)+NCOLRM);
itr = 1;
while ischar(l)
    t = regexp(l, '"(.+?)"', 'tokens');
    if length(t) < 2
        release_dates{itr} = 'NA';
        load_dates{itr} = 'NA';
    else
        release_dates{itr} = t{1}{1};
        load_dates{itr} = t{2}{1};
    end
    
    l = regexprep(l, '".+?"', '');
    
    md(itr,:) = regexp(l, ',', 'split');
    
    
    l = fgetl(fid);
    
    if mod(itr,100) == 0
        fprintf('%0.0f lines processed.\n', itr);
    end
    itr = itr + 1;
end
fclose(fid);
fprintf('%0.0f lines processed.\n', itr);

% Make the run obsnames
o = md(:,1);
md(:,1:NCOLRM) = []; % delete run ID, release date, load date.

% Add back release and load dates
md = [release_dates, load_dates, md];
vnames = [{'release_date', 'load_date'}, headers];

qt = cell2dataset(md, 'ObsNames', o, 'VarNames', vnames);
qt.spots = str2double(qt.spots);
qt.bases = str2double(qt.bases);

na = strcmpi(qt.release_date, 'NA');
qt.release_date(na) = {'Jan 00, 0000'};
qt.release_date_num = datenum(qt.release_date);
qt.release_date(na) = {'NA'};



end


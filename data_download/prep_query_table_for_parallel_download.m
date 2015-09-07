function prep_query_table_for_parallel_download( qtfile, varargin )
% Takes the query table in qtfile, filters it, and splits it into
% subtables.
% 
% Use create_srafish_calls_on_subtables to generate commands for perfoming
% parallel downloads on killdevil.

% Filter params
param.minreads = setParam(varargin, 'minreads', 4e6);
param.platform = setParam(varargin, 'platform', 'illumina');
param.processed_list = setParam(varargin, 'processed_list', []);

% Table splitting parameters
param.nsplit = setParam(varargin, 'nsplit', 50); % Number of entries per sub query table.

% Command generation parameters
param.organism = setParam(varargin, 'organism', 'Athaliana');


% Filter the query table
qt = filter_query_table(qtfile, 'minreads', param.minreads, 'platform', param.platform, ...
    'processed_list', param.processed_list);

qtfilt = [qtfile, '.filtered'];

% Create a subdirectory to work in, move the filtered query table into this
% directory, and generate all subtables in this directory.
qtdir = [regexprep(qtfilt, '\..+', ''), '_split_tables'];
mkdir(qtdir);
%system(sprintf('mv %s %s', qtfilt, qtdir));

% Split the filtered query table.
split_query_table(['../', qtfilt], qtdir, param.nsplit);

system(['tar -zcvf ', qtdir, '.tgz ', qtdir]);

% Create the download commands
%cmds = create_srafish_calls_on_subtables(qtdir, param.organism, ...
%    'donotprocess', param.processed_list);








end


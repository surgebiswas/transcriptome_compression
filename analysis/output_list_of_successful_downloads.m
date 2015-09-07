function output_list_of_successful_downloads( data_file, out_file )

load(data_file); % Assume there is a structure 's' that contains the TPM, sample IDs, and transcript IDs.

ids = s.ids;
if ~iscolumn(ids)
    ids = ids';
end

dlmcell(out_file, ids);


end


function [ s ] = read_sailfish_output( sf_file, varargin )

% If quickread is true, the TPM column is cut out, and then read in with
% matlab's 'load'.
% Note that s.table will not be a dataset but rather just a vector of TPM
% values.
quickread = setParam(varargin, 'quickread', false); 
TPMCOL = 3;

fid = fopen(sf_file);
l = fgetl(fid);

TAILRMHEADER = false;
if strcmpi(l(1), '#')
    TAILRMHEADER = false;
    s.version = l;
    s.kmer_length = fgetl(fid);
    fgetl(fid);
    s.cmd = fgetl(fid);
else
    s.version = 'NA';
    s.kmer_length = 'NA';
    s.cmd = 'NA';
end
fclose(fid);
    
if TAILRMHEADER
    system(sprintf('tail -n +6 %s > tmp.txt', sf_file));
   	fname = 'tmp.txt';
else
    fname = sf_file;
end

if quickread
    system(sprintf('cut -f %0.0f %s > tmp.txt', TPMCOL, sf_file));
    s.table = load('tmp.txt');
else
    d = dataset('file', fname, 'ReadVarNames', false, 'ReadObsNames', true);
    s.table = set(d, 'VarNames', {'Length', 'TPM', 'RPKM', 'KPKM', 'EstimatedNumKmers', 'EstimatedNumReads'});
end

if TAILRMHEADER || quickread
    !rm tmp.txt
end



end


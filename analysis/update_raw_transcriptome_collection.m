function sjoin = update_raw_transcriptome_collection( s, snew )
% Concatenates outputs of assemble_srafish_output.m script.
% s = old raw transcriptome collection struct
% snew = new raw transcriptome collection struct.


% Make sure transcript ID's are consistent between the two collections.
assert(length(s.transcript_id) == length(snew.transcript_id));

% Re sort by transcript IDs.
[s.transcript_id, sidx] = sort(s.transcript_id);
s.tpm = s.tpm(sidx,:);

[snew.transcript_id, sidx2] = sort(snew.transcript_id);
snew.tpm = snew.tpm(sidx2,:);

for i = 1 : length(s.transcript_id)
    assert(strcmpi(s.transcript_id{i}, snew.transcript_id{i}));
end

% Join the collections
sjoin.ids = [s.ids, snew.ids];
sjoin.transcript_id = s.transcript_id;

sjoin.mapped_ratio = [s.mapped_ratio, snew.mapped_ratio];
sjoin.depth = [s.depth, snew.depth];

sjoin.tpm = [s.tpm, snew.tpm];

end


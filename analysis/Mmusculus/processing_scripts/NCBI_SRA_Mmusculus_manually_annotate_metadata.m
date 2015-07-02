function NCBI_SRA_Mmusculus_manually_annotate_metadata( lY, qt, mdtablefile )
% This script will sequentially present a submissions worth of samples to
% the user.
%
% Specifically it will plot a submissions worth of samples onto the PCA. The 
% selected submission will be the submission with the largest number of
% samples that has not already been annotated in mdtablefile (an excel
% spreadsheet containing metadata). It will then be up to the user as to 
% how to annotate the data in the mdtablefile.
%
% lY = [samples x genes] expression table, log(TPM + 0.1);
% qt = query table dataset
% mdtablefile = excel spreadsheet. On the first sheet, the first column
% should be sample_id, the second column should be submission_id, and the
% remaining columns can be whatever.

load NCBI_SRA_Mmusculus_PCA_pexp_vs_eigengene_params.mat
s = lY*coef(:,1:2);

md = dataset('XLSfile', mdtablefile, 'ReadObsNames', false, 'ReadVarNames', true);

% Remove submissions in the metadata table file
usubs = unique(qt.Submission);
if isempty(unique(md.submission_id))
    qtk = qt;
else
    usubs_diff = setdiff(usubs, unique(md.submission_id));
    qtk = qt(steq(qt.Submission, usubs_idff),:);
end

% Tabulate the remaining submissions and sort by decreasing size.
tab = tabulate(qtk.Submission);
[~, sidx] = sort(cell2mat(tab(:,2)), 'descend');
tabs = tab(sidx,:);


plot(s(:,1), s(:,2), '.k');
hold on
for i = 1 : size(tabs,1)
    k = strcmpi(qt.Submission, tabs{i,1});
    qtkk = qtk(k,:);
    
    fprintf('\n\n\nSubmission: %s\n', tabs{i,1});
    plot(s(k,1), s(k,2), '.', 'Color', rand(1,3)); % enter this again if the color is not great.
    
    char(get(qtkk, 'ObsNames'))
    
    keyboard;
end




end


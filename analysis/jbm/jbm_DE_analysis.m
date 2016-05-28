function [pvals_ja, pvals_sa] = jbm_DE_analysis( y, xd )
% Compares Col-0 Mock v MeJA or BTHtreated

rep_to_use = 1;
mask_ja = strcmpi(xd.geno, 'Col-0') & xd.rep == rep_to_use & strcmpi(xd.treatment, 'MeJA');
mask_bth = strcmpi(xd.geno, 'Col-0') & xd.rep == rep_to_use & strcmpi(xd.treatment, 'BTH');
mask_mock = strcmpi(xd.geno, 'Col-0') & xd.rep == rep_to_use & strcmpi(xd.treatment, 'Mock');


pvals_ja = zeros(1, size(y,2));
pvals_sa = zeros(1, size(y,2));
for i = 1 : size(y,2)
    [~,pvals_ja(i)] = ttest2(y(mask_ja,i), y(mask_mock,i));
    [~,pvals_sa(i)] = ttest2(y(mask_bth,i), y(mask_mock,i));
end


end


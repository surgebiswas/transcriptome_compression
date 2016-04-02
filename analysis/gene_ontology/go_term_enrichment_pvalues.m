function [ pvalues ] = go_term_enrichment_pvalues( shortlistcounts, totalcounts)
% shortlistcounts = num GO terms x 1 vector of GO-term counts for the short list (interesting
%   set) of genes.
% totalcounts = '' total list of genes.


pvalues = hygecdf(shortlistcounts,max(totalcounts),...
                  max(shortlistcounts),totalcounts, 'upper');
pvalues(totalcounts <= 5) = 2;




              
end


function [ yhat ] = compression_ensemble_in_sample_predict( y, anchor_stats, cs_stats, cidx )
% For in sample prediction the generating cluster is known and is specified
% in cidx.

b_anchor = pinv(y(:,anchor_stats.S))*y;
yhat = y(:,anchor_stats.S)*b_anchor;
R = y - yhat;

for i = 1 : max(cidx)
    k = cidx == i;
    
    b_cs = pinv(R(k,cs_stats{i}.S))*R(k,:);
    yhat(k,:) = yhat(k,:) + R(k,cs_stats{i}.S)*b_cs;
end



end


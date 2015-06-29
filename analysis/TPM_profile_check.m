function [Yp, p, flagged] = TPM_profile_check( Y, varargin )
% Y = [genes x samples] TPM table.
%
% Yp = [percentiles x samples] table.

zeroprop = setParam(varargin, 'zeropropflag', 0.4);

p = 1:100;
Yp = zeros(length(p),size(Y,2));
for i = 1 : size(Y,2);
    Yp(:,i) = prctile(Y(:,i), p);
end


% Flag those with an abnormal number of zeros
flagged = mean(Yp == 0) > zeroprop;




end


function [R,d] = adaptiveqr_factor(R, d, varargin)
     [m,~] = size(R);
     if nargin > 2
       m0 = varargin{1};
     else
       [R,d] = finalize_row(R, d, 1);
       m0 = 2; 
     end
     
     % If starting row in greater than 2, perform facotrization related to the
     % starting row and below done in the iterations 2 to the starting row index
     for k = 2:(m0-1)
       [R,d] = eliminate_row(R, d, k, m0, m);
     end
     % Perform the factorization from the starting row
     for k = m0:m
       [R,d] = eliminate_row(R, d, k, k, m);
       [R,d] = finalize_row(R, d, k);
     end
end

function [R,d] = eliminate_row(R, d, k, col, m)
       Gr = R(k-1,col:m);
       d(col:m,1) = d(col:m,1) - (d(k-1).*Gr');
       R(col:m,col:m) = R(col:m,col:m) - triu(Gr'*Gr);
end

function [R,d] = finalize_row(R, d, k)
       w = sqrt(R(k,k));
       R(k,:) = R(k,:)./w;
       d(k) = d(k)/w;
end
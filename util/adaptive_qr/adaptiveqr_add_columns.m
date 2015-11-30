function [G0,T0,A] = adaptiveqr_add_columns(G0,T0,A,b,cols)
     [~,m] = size(A);
     [~,p] = size(cols);
     G0(:,m+1:m+p) = A'*cols;
     G0(m+1:m+p,m+1:m+p) = triu(cols'*cols);
     T0(m+1:m+p) = cols'*b;
     A(:,m+1:m+p) = cols;
end
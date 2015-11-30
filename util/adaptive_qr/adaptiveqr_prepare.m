function [G0,T0] = adaptiveqr_prepare(A, b)
     G0 = triu(A'*A);
     T0 = A'*b;
end
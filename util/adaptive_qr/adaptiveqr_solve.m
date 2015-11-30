function [x, R, d] = adaptiveqr_solve(G0, T0)
     [R,d] = adaptiveqr_factor(G0,T0);
     x = R\d;
end
function [ b ] = estimate_effects( x, y )

% Simple linear regression for now
% include an intercept term. if you want one.
b = pinv(x)*y; 

end


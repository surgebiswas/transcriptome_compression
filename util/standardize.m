function [ xs, mu, sig] = standardize( x, varargin)
% z-score transformation of x colum wise

centeronly = setParam(varargin, 'centeronly', false);

xs = x;
mu = zeros(1, size(xs,2));
sig = mu;
for i = 1 : size(xs,2)
    mu(i) = mean(xs(:,i));
    sig(i) = std(xs(:,i));
    
    if sig(i) == 0
        sig(i) = 1e-10;
    end
    
    if centeronly
        xs(:,i) = xs(:,i) - mean(xs(:,i));
    else
        xs(:,i) = (xs(:,i) - mean(xs(:,i)))/sig(i);
    end
    
end

end


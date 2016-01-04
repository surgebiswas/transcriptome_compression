function w = gaussian_kernel_diagonal( x, x0, lambda )

% lambda should be a [1 x size(x,2)] vector representing the diagonal
% entries of the diagonal [size(x,2) x size(x,2)] covariance matrix.
% if lambda is a scalar then it is assumed that the covariance matrix is
% given by lambda*eye(size(x,2))

c = bsxfun(@minus, x, x0).^2;

if size(x,2) > 1 && length(lambda) == 1
    w = exp( -sum(c,2)/(2*lambda) );
else
    w = exp(-sum(bsxfun(@rdivide,c,2*lambda),2));
end






end


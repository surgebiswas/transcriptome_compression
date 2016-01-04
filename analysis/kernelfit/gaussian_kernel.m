function [ w ] = gaussian_kernel( x, x0, lambda )


c = bsxfun(@minus, x, x0);

w = zeros(size(x,1),1);

for i = 1 : length(w)
    w(i) = exp( -c(i,:)*(lambda\c(i,:)')/2 );%/sqrt( ((2*pi)^size(x0,2))*det(lambda) );
end



end


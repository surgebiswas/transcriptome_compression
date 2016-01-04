function [ xu ] = unstandardize( x, mu, sig )

xu = bsxfun(@plus, bsxfun(@times, x, sig), mu);


end


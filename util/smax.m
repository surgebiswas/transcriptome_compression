function [ sm ] = smax( x )

sm = exp(x-logsumexp(x));

end


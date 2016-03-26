function er = erank( x )

S = svd(x);

p = S/sum(abs(S));

er = exp( -sum( p.*log(p) ) );


end


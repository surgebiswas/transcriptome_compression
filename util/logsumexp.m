function s=logsumexp(x)
y=max(x);
s=y+log(sum(exp(x-y)));

end


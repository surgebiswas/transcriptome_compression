function invwish_optim_delta(n, d )
% Empirical bayes (type II max likelihood) estimate of the d.o.f of a
% inverse wishart.

loglike = @(x) (x*d/2)*(log(x+n) - log(x)) + (n*d/2)*log(x + n) + logMvGamma(x/2, d) - logMvGamma((x+n)/2,d); 

q = linspace(d-1+0.1,1000*d,100);
fv = zeros(size(q));
gfv = fv;
for k = 1 : length(fv)
fv(k) = loglike(q(k));
gfv(k) = grd(q(k),n,d);
end
subplot(1,2,1)
plot(q,fv, '-k')
subplot(1,2,2)
plot(q,gfv, '-k')


    function g = grd(x,n,d)
        t1 = (d/2)*(log(x + n) - log(x));
        t2 = 0;
        for i = 1 : n
            for j = 1 : n
                t2 = t2 + 1/(x + n - i - 1 - j);
            end
        end
        g = t1 - t2;
    end

end


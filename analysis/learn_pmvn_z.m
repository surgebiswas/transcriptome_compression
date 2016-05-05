function [ z ] = learn_pmvn_z( t, mu, Sigma, offset )

% t is observations x variables
% mu is 1 x variables
% sigma is variables x variables covariance matrix
% offset is observations x 1
%
% Model: z ~ Normal(mu, Sigma); t ~ Poisson(exp(z)*offset)
%offset = repmat(offset,1,size(t,2));
z = log(t+1) - repmat(log(offset),1, size(t,2)); %mvnrnd(mu, Sigma, size(t,1));

%optim_opts = optimoptions('fminunc','Algorithm','trust-region', ...
%    'SpecifyObjectiveGradient',true, 'HessianFcn', 'objective');

delta = Inf;
o = -Inf;
while delta > 1e-6
    % update a column of z at a time.
    for j = 1:size(z,2)
        [mu_cond, sigma_cond] = get_conditional_normal_params(z, mu, Sigma, j);
        
%         ofun = @(x) zfun(x, t(:,j), mu_cond, sigma_cond, offset);
%         zup = fminunc(ofun, z(:,j), optim_opts);
%         z(:,j) = zup;

        [f,g,H] = zfun(z(:,j), t(:,j), mu_cond, sigma_cond, offset);
        z(:,j)  = z(:,j) - g./H; % newton raphson update
        
        
        % Using Lambert W function causes overflow issues. DEPRECATED
        % z(:,j) = mu_cond + sigma_cond*t(:,j) - ...
        %    lambertw(0, sigma_cond*offset.*exp(mu_cond + sigma_cond*t(:,j)) );
        
    end
    old_o = o;
    o = objective(t, z, offset, mu, Sigma);
    delta = o - old_o;
    fprintf('Objective: %0.6f\n', o);
end

    function o = objective(t, z, o, mu, S)
        lambda = exp(z).*repmat(o, 1, size(z,2));
        
        poiss_comp = sum(sum( t.*log(lambda) - lambda  ));
        gauss_comp = sum(logmvnpdf(z, mu, S));
        o = poiss_comp + gauss_comp;
    end


    function [mc,sc] = get_conditional_normal_params(z, m, s, i)
        ci = [1:i-1,i+1:size(m,2)];
        sub = bsxfun(@minus, z(:,ci), m(ci));
        
        cinv = s(ci,ci)\s(ci,i);
        
        mc = m(i) + sub*cinv;
        sc = s(i,i) - s(i,ci)*cinv;
    end

    function [f,g,H] = zfun(zk, tk, mc, sc, o)
        a = exp(zk).*o;
        zc = zk - mc;
        sci = 1/sc;
        
        f = tk.*zk - a - 0.5*sci*zc.^2;
        
        if nargout > 1
            g = tk - a - zc*sci;
            
            if nargout > 2
                H = -a - sci;
                %H = spdiags(h, 0, length(h), length(h));
            end
        end
        
    end




end


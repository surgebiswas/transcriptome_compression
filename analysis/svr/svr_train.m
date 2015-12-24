function [ model, perfstats ] = svr_train( x, y, varargin )
% svr_train - Training for support vector regression.

nsearch = setParam(varargin, 'nsearch', 10);
nsearch_fine = setParam(varargin, 'nsearch_fine', 10);
cvi = setParam(varargin, 'crossvalind', crossvalind('Kfold', size(x,1), nsearch));
C = setParam(varargin, 'C', logspace( log10(2*3*std(y)) - 2, log10(2*3*std(y)) + 2, nsearch));
gamma = setParam(varargin, 'gamma', logspace(-6, 3, 10));
epsilon = setParam(varargin, 'epsilon', set_default_epsilon(x,y,nsearch) );
finesearch = setParam(varargin, 'finesearch', true);


[tp, tps] = cross_validate(x, y, C, gamma, epsilon, cvi);

if finesearch
    [c_fine, gamma_fine, eps_fine] = build_fine_search_schedule(tp, tps, C, gamma, epsilon, nsearch_fine);
    [tp_fine, tps_fine] = cross_validate(x, y, c_fine, gamma_fine, eps_fine, cvi);
    
    [ci,gi,ei] = find_optimal_parameter(tp_fine,tps_fine);
else
    [ci,gi,ei] = find_optimal_parameter(tp,tps);
end


% Train on full data.
model = svmtrain(y, x, sprintf('-s 3 -t 2 -p %0.20f -g %0.20f -c %0.20f -q', epsilon(ei), gamma(gi), C(ci)));
perfstats.ci = ci;
perfstats.gi = gi;
perfstats.ei = ei;
perfstats.test_perf = tp;
perfstats.test_perf_se = tps;

if finesearch
    perfstats.test_perf_fine = tp_fine;
    perfstats.test_perf_se_fine = tps_fine;
end



    function [cf,gf,ef] = build_fine_search_schedule(tp,tps,c,g,e,nf)
        [cidx,gidx,eidx] = find_optimal_parameter(tp,tps);
        
        
        cf = set_schedule(c,cidx,nf);
        gf = set_schedule(g,gidx,nf);
        ef = set_schedule(e,eidx,nf);
        
    end

    function s = set_schedule(p, pidx, nf)
        
        if pidx == 1
            pmin = p(1)/10;
            pmax = p(2);
        elseif pidx == length(p)
            pmin = p(end-1);
            pmax = p(end)*10;
        else
            pmin = p(pidx-1);
            pmax = p(pidx+1);
        end
        
        s = logspace(log10(pmin),log10(pmax), nf);
        
    end

    function eps = set_default_epsilon(x,y,n)
        yhat = mean(y(knnsearch(x,x, 'K', round(0.01*length(y)))),2);
        sigma = std( (y - yhat) );
        
        eps = logspace( log10(sigma) - 2, log10(sigma) + 2, n);
    end

    function [test_perf, test_perf_std] = cross_validate(x, y, C, gamma, epsilon, cvi)
        test_perf = zeros(length(C), length(gamma), length(epsilon));
        test_perf_std = test_perf;
        for c = 1 : length(C)
            for g = 1 : length(gamma)
                for e = 1 : length(epsilon)
                
                    %%%% INNER CROSS VALIDATION LOOP %%%%
                    cvperf = zeros(max(cvi),1);
                    for i = 1 : max(cvi)
                        xtrain = x(cvi~=i,:);
                        ytrain = y(cvi~=i);

                        xtest = x(cvi==i,:);
                        ytest = y(cvi==i);

                        m = svmtrain(ytrain, xtrain, sprintf('-s 3 -t 2 -p %0.20f -g %0.20f -c %0.20f -q', epsilon(e), gamma(g), C(c)));

                        [yhat,acc,~] = svmpredict(ytest, xtest, m, '-q');
                        cvperf(i) = acc(2); 
                    end
                    test_perf(c,g,e) = mean(cvperf); % Store average MSE across all folds.
                    test_perf_std(c,g,e) = std(cvperf)/sqrt(length(cvperf)); % Standard error of mse's.
                    %%%% INNER CROSS VALIDATION LOOP %%%%
                    
                    
                end
                
            end
            
            if true; fprintf('Cross validation %0.2f%% complete.\n', 100*c/length(C)); end
        end
    end
    
    function [c_star_idx, gamma_star_idx, epsilon_star_idx] = find_optimal_parameter(test_performance, test_performance_std)
        % In this current implementation we will just look what gives the
        % best test_performance. We assume test_performance is given by
        % mean squared error (lower is better).
        % 
        % test_performance and test_performance_std are |C| x |gamma| x |epsilon| 
        % matrices representing test set performances over schedules of C, 
        % gamma and epsilon respectively.
        
        bc = zeros(size(test_performance,3),1);
        bg = zeros(size(test_performance,3),1);
        vals = bc;
        for i = 1 : size(test_performance,3)
            [vals(i),bc(i),bg(i)] = max2d(-test_performance(:,:,i));
        end
        
        [~,epsilon_star_idx] = max(vals);
        c_star_idx = bc(epsilon_star_idx);
        gamma_star_idx = bg(epsilon_star_idx);
        
    end
    
    



end


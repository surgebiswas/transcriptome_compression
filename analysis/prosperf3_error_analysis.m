function prosperf3_error_analysis( results, sets, tids, lY, qt )

tp = results.target_proc;
pp = results.pred_proc;
tg = results.target_gene;
pg = results.pred_gene;
cvi = results.cvi;

% Build geneset membership matrix
disp('Building tr. program membership matrix.');
Q = false(size(tg,2), length(sets));
for i = 1 : length(sets)
    Q(:,i) = steq(tids, sets{i});
end



% Programs first.
if false
    disp('Analyzing error for tr. programs');
    tpa = subadjust(tp, qt.Submission); 
    ppa = subadjust(pp, qt.Submission);

    % remove low number programs.
    pn = sum(Q);
    tpa(:,pn<10) = [];
    ppa(:,pn<10) = [];
    Q(:,pn<10) = [];


    pcc = rsq_and_slope(tpa, ppa);
    uvar = var(tpa - ppa)./var(tpa);


    vexp = var(tpa); % Variance of program
    ab = zeros(1, size(tpa,2));
    for j = 1 : size(tpa,2)
        ab(j) = mean(mean(lY(:,Q(:,j)))); % abundance of program
    end
    psize = sum(Q);


    % vexp(:, pn<10) = [];
    % ab(:, pn<10) = [];
    % psize(:,pn<10) = [];

    xlabs = {'log(Expression variance)', 'Avg. member abundance', 'log(Program size)'};
    xvals = {log(vexp), (ab), log(psize)};

    ylabs = {'PCC', 'Unexplained variance'};
    yvals = {pcc, uvar};


    idx = 1;
    for i = 1 : length(yvals)
        for j = 1 : length(xvals)
            subplot(2,3,idx);
            if j == 1
                if i == 2; 
                    err_analysis_plot(xvals{j}, yvals{i}, xlabs{j}, ylabs{i})
                else
                    err_analysis_plot(xvals{j}, yvals{i}, [], ylabs{i})
                end
            else
                if i == 2; 
                    err_analysis_plot(xvals{j}, yvals{i}, xlabs{j}, [])
                else
                    err_analysis_plot(xvals{j}, yvals{i}, [], [])
                end
            end
            idx = idx + 1;
        end
    end

    plotSave('figures/error_analysis/error_v_variables.png');


    xx = [log(vexp'), ab', log(psize')];
    st = regstats(log(uvar'), xx, 'linear');
    phat = xx*st.beta(2:end) + st.beta(1);

    figure;
    subplot(2,1,1);
    err_analysis_plot(pcc, uvar, 'PCC', 'Unexplained variance');


    subplot(2,1,2);
    err_analysis_plot(phat', log(uvar), 'Predicted log(Unexp. Var.)', 'Actual log(Unexp. Var.)');
    plotSave('figures/error_analysis/error_prediction.png');
    close all
end

% Genes
if true
    disp('Analyzing error for genes');
    tga = subadjust(tg, qt.Submission); 
    pga = subadjust(pg, qt.Submission);
    
    pcc = rsq_and_slope(tga, pga);
    uvar = var(tga - pga)./var(tga);
    
    vexp = var(tga);
    nprog = sum(Q,2)';
    ab = mean(lY);
    
    xlabs = {'log(Expression variance)', 'Avg. abundance', 'Num. programs'};
    xvals = {log(vexp), (ab), nprog + 0.5*rand(size(nprog)) };

    ylabs = {'PCC', 'Unexplained variance'};
    yvals = {pcc, uvar};
    
    idx = 1;
    for i = 1 : length(yvals)
        for j = 1 : length(xvals)
            subplot(2,3,idx);
            if j == 1
                if i == 2; 
                    err_analysis_plot(xvals{j}, yvals{i}, xlabs{j}, ylabs{i}, 'msize', 1)
                    %v = axis;
                    %axis([v(1) v(2) 0 1]);
                else
                    err_analysis_plot(xvals{j}, yvals{i}, [], ylabs{i}, 'msize', 1)
                end
            else
                if i == 2; 
                    err_analysis_plot(xvals{j}, yvals{i}, xlabs{j}, [], 'msize', 1)
                    %v = axis;
                    %axis([v(1) v(2) 0 1]);
                else
                    err_analysis_plot(xvals{j}, yvals{i}, [], [], 'msize', 1)
                end
            end
            
            if j == 3;
                set(gca, 'XTick', 0:2:max(xvals{3}));
                v = axis;
                axis([-0.5 v(2:end)]);
            end
            
            idx = idx + 1;
        end
    end
    plotSave('figures/error_analysis/error_v_variables_genes.png');
    close
    
    xx = [log(vexp'), ab', nprog'];
    st = regstats(log(uvar'), xx, 'linear');
    phat = xx*st.beta(2:end) + st.beta(1);

    figure;
    subplot(2,1,1);
    err_analysis_plot(pcc, uvar, 'PCC', 'Unexplained variance', 'msize', 1);

    subplot(2,1,2);
    err_analysis_plot(phat', log(uvar), 'Predicted log(Unexp. Var.)', 'Actual log(Unexp. Var.)', 'msize', 1);
    plotSave('figures/error_analysis/error_prediction.png');
    close all

end




% mean abundance
% pcc = zeros(max(cvi), size(tp,2));
% rmse = pcc;
% vartrain = pcc;
% vartest = pcc;
% ab = pcc;
% psize = sum(Q);
% for i = 1 : max(cvi)
%     
%     vartrain(i,:) = var(tp(cvi ~= i,:));
%     vartest(i,:) = var(tp(cvi == i,:));
%     
%     
%     for j = 1 : size(tp,2)
%         ab(i,j) = mean(mean(lY(cvi~=i,Q(:,j))));
%     end
%     
%     
%     pcc(i,:) = rsq_and_slope(tpa(cvi==i,:), ppa(cvi==i,:));
%     rmse(i,:) = var( tpa(cvi==i,:) - ppa(cvi==i,:) )./var( tpa(cvi==i,:) );
%     
%     
%     disp(i);
% end








end


function analyze_posterior_predictive_dist( Y, qt, tids, sets )

path(genpath('~/GitHub/tradict'), path)
[ytrain, ytest, ktrain] = partition_data(Y', qt, 0.1);
nmax = 20;

T_train = ytrain.*repmat(qt.spots(ktrain)/1000000,1, size(ytrain,2));
o_train = qt.spots(ktrain)/1000000;

T_test = ytest.*repmat(qt.spots(~ktrain)/1000000,1, size(ytest,2));
o_test = qt.spots(~ktrain)/1000000;


model = tradict_train( T_train, o_train, tids, sets );

rng(23);
idx = randsample(size(T_test,1), nmax);
T_m = T_test(:,model.S);
pred = tradict_predict( T_m(idx,:), o_test(idx,:), model, 'sample_posterior', true );





% Actual values
z = lag_dataset(T_test, o_test, 'priors', model.lag_priors);
zs = standardize(z, 'mu', model.train_mu, 'std', model.train_sig);
s = zs*model.geneset.coef;

zs = z(idx,:);
ss = s(idx,:);


[ cs, ac ] = compare_post_pred_dist( ss, pred.programs.post_smpl );

gidx = randsample(size(zs,2), 1000);
ps = pred.genes.post_smpl;
for i =1 : length(ps)
    ps{i} = ps{i}(:,gidx);
end
[ csg, acg ] = compare_post_pred_dist( zs(:,gidx), ps );


figure;
subplot(1,2,2)
plot([0 1], [0 1], '--k', 'LineWidth', 2);
hold on
serr = sqrt(0.01*cs.*(1 - cs));
plot(cs, ac, '-r', 'LineWidth', 2);
jbfill(cs, ac+serr, ac-serr, 'r', 'r', true, 0.3);
xlabel('Credible interval size', 'FontSize', 18);
set(gca, 'FontSIze', 14);
title('tr. programs', 'FontSize', 20);
set(gca, 'XTick', 0:0.2:1);
set(gca, 'XTickLabel', 100*[0:0.2:1]);
set(gca, 'YTickLabel', 100*[0:0.2:1]);
axis square

subplot(1,2,1);
plot([0 1], [0 1], '--k', 'LineWidth', 2);
hold on
serr = sqrt(0.01*cs.*(1 - cs));
plot(csg, acg, '-g', 'LineWidth', 2);
jbfill(csg, acg+serr, acg-serr, 'g', 'g', true, 0.3);
xlabel('Credible interval size', 'FontSize', 18);
set(gca, 'FontSIze', 14);
ylabel('Perc. of test set contained', 'FontSize', 18);
title('genes', 'FontSize', 20);
set(gca, 'XTick', 0:0.2:1);
set(gca, 'XTickLabel', 100*[0:0.2:1]);
set(gca, 'YTickLabel', 100*[0:0.2:1]);
axis square


plotSave('figures/credible_interval_validation.png');
close




% ss_true_std = zeros(size(ss));
% ss_dist_std = cell(size(ss));
% xi = linspace(-4,4,100);
% for i = 1 : length(idx)
%     for j = 1 : size(ss,2)
% %         mm = mean(pred.programs.post_smpl{i}(:,j));
% %         st = std(pred.programs.post_smpl{i}(:,j));
% %         
% %         
% %         x = (pred.programs.post_smpl{i}(:,j) - mm)/st;
% %         ss_dist_std{i,j} = ksdensity(x, xi);
% %         ss_dist_std{i,j} = ss_dist_std{i,j}/sum(ss_dist_std{i,j});
% %         ss_true_std(i,j) = (ss(i,j) - mm)/st;
%     end
% end



% 
% figure
% nn = hist(ss_std(:), 100);
% mn = max(nn);
% hold on
% F = zeros(length(idx)*size(ss,2), length(xi));
% ii = 1;
% for i = 1 : length(idx)
%     for j = 1 : size(ss,2)
%         F(ii,:) = ss_dist_std{i,j}*mn./max(ss_dist_std{i,j});
%         ii = ii + 1;
%         %plot(xi, ss_dist_std{i,j}*mn./max(ss_dist_std{i,j}), '-r', 'LineWidth', 0.5)
%     end
% end



end


clear
d = dataset('file', '~/GitHub/tradict/ncbi_sra_growth.txt', 'ReadObsNames', false, 'ReadVarNames', true);
d(end,:) = [];
mindate = min(datenum(d.DateAsNumber));
maxdate = max(datenum(d.DateAsNumber));
x = datenum(d.DateAsNumber); x= x- mindate;
y = cumsum(d.count);

[ux,ia] = unique(flipud(x));
yud = flipud(y);
yudu = yud(ia);

plot(ux, yudu, 'ok')

hold on 
F = @(x,xdata)x(1)*exp(-x(2)*xdata) + x(3)*exp(-x(4)*xdata);
x0 = [1 1 1 0];
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,ux,yudu);
plot(ux, F(x,ux), '-b', 'LineWidth', 2)


expVal = 185;
query = [max(ux):10:(max(ux)+expVal)];
plot(query, F(x,query), '--b', 'LineWidth', 2)
plot(query(end), F(x,query(end)), 'ob', 'LineWidth', 3)

axis square;
axis tight;
set(gca, 'XTick', fliplr((max(ux)+expVal-10):-365:min(ux)))
set(gca, 'XTickLabel', flipud(cellstr(datestr((maxdate+185):-365:mindate))));
ytick = round(F(x,fliplr((max(ux)+expVal-10):-365:min(ux))'));
set(gca, 'YTick', ytick);
ylabel('Number of {\itA. thaliana } RNA-Seq datasets in NCBI SRA')
grid on
rotateXLabels(gca, 36);
set(gca, 'FontSize', 8)
plotSave('/Users/sbiswas/Documents/matlab/src/rampseq/figures/At_NCBI_SRA_growth.png')
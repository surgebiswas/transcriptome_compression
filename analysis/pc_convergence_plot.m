function  pc_convergence_plot( s, qt, pexp )


dn = qt.release_date_num;
udn = unique(qt.release_date_num);


res = 1000;

% Generate density estimates of the full data.
fins = zeros(res, size(s,2));
grids = zeros(res, size(s,2));
for i = 1 : size(s,2)
    grids(:,i) = linspace(min(s(:,i)), max(s(:,i)), res);
    fins(:,i) = ksdensity(s(:,i), grids(:,i));
end

if false
    R = zeros(size(s,2), length(udn));
    parfor (i = 1 : length(udn),  6)
        mask = dn <= udn(i);
        
        r = zeros(size(s,2),1);
        for j = 1 : size(s,2)
            dd = ksdensity(s(mask,j), grids(:,j));
            r(j) = corr(dd, fins(:,j));
        end
        R(:,i) = r;
        disp(i);
    end

    save('NCBI_SRA_PC_convergence_plot.mat', 'R');
else 
    load('NCBI_SRA_PC_convergence_plot.mat');
end

% Plot
figure;
hold on
cm = flipud(repmat(linspace(0,0.95,res)', 1, 3));
colormap(cm);
for i = size(s,2) : -1 : 1
    cm_idx = round(size(cm,1)*pexp(i)/max(pexp));
    
    plot(udn, R(i,:), '-', 'Color', cm(cm_idx,:), 'LineWidth', 2)
end
axis square
axis tight
set(gca,'Layer','top')
v = axis;

xtick = fliplr(max(udn) : -365 : min(udn)); %max(udn)-10*365 : 365 : max(udn);
set(gca, 'XTick', xtick);
set(gca, 'XTickLabel', datestr(xtick, 'mmm-yyyy'));
axis(v);
box on;
set(gca, 'FontSize', sf.axis_tick_labels-2);


try; rotateXLabels(gca, 20); end
ylabel('Pearson correlation', 'FontSize', sf.axis_labels);

c  = colorbar;
caxis([pexp(size(s,2)) max(pexp)])
plotSave('figures/pca/pc_convergence.png');
close




end


clear;
rng('default')
cd ~/GitHub/data/transcriptome_compression/fig1/

path(genpath('~/GitHub/transcriptome_compression'), path)
bmus = [-2 0; 1 -1; 0 1; -3 2; 1 0];

b = [];
for i = 1 : size(bmus,1)
    b = [b; randn(5,2)*0.5 + repmat(bmus(i,:), 5, 1)];
end
b = b';

markers = randn(50,2);

y = standardize(markers*b);
ygs = standardize(markers*bmus');

ri=hclust(y);
y = y(ri,:);
ygs = ygs(ri,:);
markers = markers(ri,:); %markers(:,1) = -markers(:,1);


% Full gene expression matrix
save_fig1_subheatmap(y, 'gene', 'full_gene_heatmap.png');

% Gene set matrix
save_fig1_subheatmap(ygs, 'geneset', 'geneset_heatmap.png');


% Separated 
mreps = [5, 14]; 
save_fig1_subheatmap(y(:,1:mreps(1)-1), 'gene', 'sub_gene_heatmap_1.png');
save_fig1_subheatmap(y(:,mreps(1)+1:mreps(2)-1), 'gene', 'sub_gene_heatmap_2.png');
save_fig1_subheatmap(y(:,mreps(2)+1:end), 'gene', 'sub_gene_heatmap_3.png');

save_fig1_subheatmap(y(:,mreps(1)), 'gene', 'marker1_heatmap.png');
save_fig1_subheatmap(y(:,mreps(2)), 'gene', 'marker2_heatmap.png');

% Measurements -- decoding
subidx = 1:2:size(y,1);
save_fig1_subheatmap(-markers(subidx,1), 'gene', 'measured_marker1.png'); % negative because we're using col 5 in the full exp mat.
save_fig1_subheatmap(markers(subidx,2), 'gene', 'measured_marker2.png');


ygsh = standardize(markers*bmus' + 0.3*randn(size(ygs)));
save_fig1_subheatmap(ygsh(subidx,:), 'geneset', 'geneset_predicted_heatmap.png');

yh = standardize(markers*b + 0.3*randn(size(y)));
save_fig1_subheatmap(yh(subidx,:), 'gene', 'full_gene_predicted_heatmap.png');

% Colorbars
figure
colormap(((cbrewer('seq', 'YlOrRd', 64))))
c = colorbar;
set(c, 'TickLabels', []);
set(c, 'TickLength', [0 0])
axis off
plotSave('TPM_colorbar.png');
close

figure
colormap(parula)
c = colorbar;
set(c, 'TickLabels', []);
set(c, 'TickLength', [0 0])
axis off
plotSave('geneset_colorbar.png');
close


figure
colormap(prgn)
c = colorbar;
set(c, 'TickLabels', []);
set(c, 'TickLength', [0 0])
axis off
plotSave('covar_colorbar.png');
close



% Distribution plots
if false
    figure;
    ar = [2.5,1,1];
    npdf = @(x,m,s) normpdf(x,m,s)/normpdf(m,m,s); % scaled.
    ppdf = @(x,l) ((l.^x)*exp(-l)./gamma(x + 1))/((l.^floor(l))*exp(-l)./gamma(floor(l) + 1));

    clf
    subplot(3,1,1)
    x = -4:0.01:4;
    n = npdf(x,0,1);
    rvals = [-1.4, 0.8];
    hold on
    plot([rvals(1) rvals(1)], [0 npdf(rvals(1),0,1)], '-r', 'LineWidth', 3) 
    plot([rvals(2) rvals(2)], [0 npdf(rvals(2),0,1)], '-r', 'LineWidth', 3) 
    ah = area(x,n);
    set(ah, 'FaceColor', [0.4 0.4 0.4])
    set(ah, 'FaceAlpha', 0.3)
    set(ah, 'LineWidth', 3)
    axis off
    daspect(ar)


    subplot(3,1,2); cla;  hold on
    n1 = npdf(x,rvals(1),0.6);
    n2 = npdf(x,rvals(2),0.4);
    rvals2 = rvals;
    %plot([rvals2(1) rvals2(1)], [0 npdf(rvals2(1),rvals(1), 0.6)], '-r', 'LineWidth', 3) 
    %plot([rvals2(2) rvals2(2)], [0 npdf(rvals2(2),rvals(2), 0.4)], '-r', 'LineWidth', 3) 

    ah1 = area(x,n1);
    set(ah1, 'FaceColor', [0 0 1])
    set(ah1, 'FaceAlpha', 0.3)
    set(ah1, 'LineWidth', 3)
    ah2 = area(x,n2);
    set(ah2, 'FaceColor', [1 0 0])
    set(ah2, 'FaceAlpha', 0.3)
    set(ah2, 'LineWidth', 3)
    axis off
    daspect(ar)



    subplot(3,1,3); cla; hold on

    p1 = ppdf(exp(x), exp(rvals2(1)));
    p2 = ppdf(exp(x), exp(rvals2(2)));

    ah21 = area(exp(x),p1);
    set(ah21, 'FaceColor', [0 0 1])
    set(ah21, 'FaceAlpha', 0.3)
    set(ah21, 'LineWidth', 3)
    ah22 = area(exp(x),p2);
    set(ah22, 'FaceColor', [1 0 0])
    set(ah22, 'FaceAlpha', 0.3)
    set(ah22, 'LineWidth', 3)

    axis([0 7 0 1.15]);
    axis off
    daspect([1.96 1 1])

    plotSave('fig1_distributions.png');
    close
end


% MVN illustration plot
if true
    rng(2)
    S = [1 -0.4 0.85; -0.4 1 -0.6; 0.85 -0.6 1];
    mu = zeros(1,3)+1;
    
    z = mvnrnd(mu, S, 100);
    t = poissrnd(exp(z));
    [~,sidx] = sort(t(:,1), 'descend');
    
    figure;
    for i = 1 : size(t,2)
        subplot(3,1,i)
        b = bar( (t(sidx,i)) );
        set(b, 'FaceColor', 'k')
        set(b, 'BarWidth', 1)
        axis tight;
        axis off
    end
    plotSave('count_bars.png');
    close
    
    
    
    f = figure;
    clf
    [h, ax, bax, p, pax] = plotmatrix(z, '.k');
    offset = 1;
    for i = 1 : length(h)
        set(h(i,i), 'XData', []);
        set(h(i,i), 'YData', []);
        set(pax(i), 'Visible', 'off')
        v = axis( ax(i,i) );
        axis(ax(i,i), [-5 5 v(3) v(4)]);
        set(ax(i,i), 'Visible', 'off')
        set(pax(i), 'Xlim', [min(z(:,i)) - offset, max(z(:,i)) + offset]);
        set(p(i),'EdgeColor',[0.5 0.5 0.5],'FaceColor',[0.5 0.5 0.5]);
        set(p(i), 'BinEdges', [-3:0.5:3]);
    end


    for i = 1 : size(z,2)
        for j = i + 1 : size(z,2)
            set(ax(i,j), 'Visible', 'off');
            set(h(i,j), 'Marker', 'o');
            set(h(i,j), 'MarkerSize', 8);
            set(h(i,j), 'MarkerFaceColor', 'k');
            set(ax(i,j), 'Xlim', [min(z(:,j)) - offset, max(z(:,j)) + offset]);
            set(ax(i,j), 'Ylim', [min(z(:,i)) - offset, max(z(:,i)) + offset]);
        end
    end

    for j = 1 : size(z,2)
        for i = j + 1 : size(z,2)
            set(ax(i,j), 'Visible', 'off');
            set(h(i,j), 'Xdata', []);
            set(h(i,j), 'Ydata', []);
        end
    end
    
    
    axis([ax(:)], 'square')
    axis([pax(:)], 'square')
    % Save manually
    
    
    figure
    imagesc(S, [-1 1]); colormap(prgn)
    c = colorbar;
    set(c, 'TickLength', [0 0])
    set(c, 'TickLabels', []);
    axis square
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca, 'TickLength', [0 0]);
    add_heatmap_gridlines(S); 
    set(gca, 'LineWidth', 5)
    plotSave('fig1_covar.png');
    close
    
    
    
end
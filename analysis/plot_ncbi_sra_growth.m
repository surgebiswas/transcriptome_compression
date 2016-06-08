function plot_ncbi_sra_growth( qt )
% qt = query table read in from read_ncbi_sra_query_table.

EXTRAPVAL = 0; % number of days to extrapolate out to
sf = get_standard_figure_font_sizes;

qt( strcmpi(qt.release_date, 'NA'), :) = [];


dn = qt.release_date_num;
mindn = min(dn);

if true
    dnu = unique(dn) - mindn;
    nlib = zeros(length(dnu),1);
    for i = 1 : length(nlib)
        nlib(i) = sum( dn <= dnu(i)+mindn);
    end
else
    load('tmp.mat');
end


%F = @(x,xdata)x(1)*exp(-x(2)*xdata) + x(3)*exp(-x(4)*xdata);
% F = @(x,xdata)x(1)*exp(-x(2)*xdata);
% x0 = [1 1];
% x = lsqcurvefit(F,x0,dnu,nlib);


% Exclude the first 40% of data points when estimating the growth trend.
b = estimate_growth_trend(dnu,nlib,round(length(dnu)*0.4));

plot(dnu, nlib, 'ok');
hold on;
plot(dnu, exp(dnu*b(2) + b(1)), '-b', 'LineWidth', 2);

addx = [max(dnu):10:max(dnu)+EXTRAPVAL];
plot(addx, exp(addx*b(2) + b(1)), '--b', 'LineWidth', 2);
plot( addx(end), exp(addx(end)*b(2) + b(1)), 'ob', 'LineWidth', 2);

axis square;
axis tight;
box on


xtick = (max(dnu)+EXTRAPVAL):-365:0;
xticklabel = datestr(xtick + mindn);


set(gca, 'XTick', fliplr(xtick));
set(gca, 'XTickLabel', flipud(cellstr(xticklabel)));
try; rotateXLabels(gca, 36);end
ytick =  floor(exp(xtick*b(2) + b(1)));

NYTICKSKIP = 3; ytick = [ytick(1:end-NYTICKSKIP), 0];
set(gca, 'YTick', fliplr(ytick));
set(gca, 'YTickLabel', cellstr(num2str(fliplr(ytick)'))  );
grid on;
set(gca,'GridLineStyle','-');
set(gca, 'FontSize', sf.axis_tick_labels-4);
ylabel('Num. HQ transcriptomes', 'FontSize', sf.axis_labels);

v = axis;
axis([v(1), v(2), v(3) - v(4)*0.01, v(4)]);


    function b = estimate_growth_trend(x,y,numexclude)
        % Estimates trend in the log-scale of y.
        
        if ~iscolumn(x)
            x = x';
        end
        
        if ~iscolumn(y)
            y = y';
        end
        
        xe = x; xe(1:numexclude) = [];       
        ye = y; ye(1:numexclude) = [];
        
        b = pinv([ones(size(xe,1),1), xe])*log(ye);
        
    end






end


function err_analysis_plot( x, y, xlab, ylab, varargin)

sp = setParam(varargin, 'usespline', true);
msize = setParam(varargin, 'msize', 10);
sr = setParam(varargin, 'splinerange', 'minmax');

plot(x,y, '.k', 'MarkerSize', msize);
hold on

h = mean(diff(sort(x,'ascend')));
pw = 1/(1 + (h^3)/6);

if strcmpi(sr, 'minmax')
    pp = csaps(x,y,0.4); %spline(x,y);
elseif strcmpi(sr, '3sig')
    xx = linspace(mean(x)-3*std(x), mean(x)+3*std(x), 100);
    pp = csaps(x,y,0.4, xx);
end
    
% xrange = linspace(min(x), max(x), 100);
% yf = ppval(pp, xrange);

if sp
    if strcmpi(sr, '3sig')
        plot(xx, pp, '-r', 'LineWidth', 2);
    else
        fnplt(pp,'-r');
    end
end

axis tight;
axis square
buffer_axis;

sf = get_standard_figure_font_sizes;
if ~isempty(xlab)
    %xlabel(xlab, 'FontSize', sf.axis_labels);
    xlabel(xlab);
end

if ~isempty(ylab)
    %ylabel(ylab, 'FontSize', sf.axis_labels);
    ylabel(ylab);
end

v = axis;
tx = v(1) + 0.03*(v(2) - v(1));
ty = v(3) + 0.94*(v(4) - v(3));
text(tx, ty, ['\rho = ', sprintf('%0.2f', corr(x',y', 'type', 'spearman'))]);

%set(gca, 'FontSize', sf.axis_tick_labels);



end


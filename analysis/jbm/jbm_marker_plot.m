function  jbm_marker_plot( y, yhat, xd, treatment, genos, ylab, varargin)

rep_to_use = setParam(varargin, 'rep_to_use', 1);

mask = steq(xd.geno, genos) & xd.rep == rep_to_use & (strcmpi(xd.treatment, treatment) | strcmpi(xd.treatment, '-'));

geno2color = containers.Map;
cm = lines;
geno2color('coi1-16') = cm(1,:);% Blue %'r';
geno2color('npr1-1') = cm(2,:); % red  %'b';
geno2color('Col-0') = cm(3,:); % yellow 'k';
geno2color('35S:HopBB1') = cm(4,:); % purple   'g';
geno2color('eds16-1') = cm(5,:); % green 'c';


xdm = xd(mask,:);
xdm.time = str2double(xdm.time);
ym = y(mask,:);
yhatm = yhat(mask,:);


times = [0 1 5 8];
for i = 1 : length(genos)
    
    q = strcmpi(xdm.geno, genos{i});
    
    pv = [];
    pvh = [];
    for j = 1 : length(times)
        pv = [pv, mean(ym(q&xdm.time == times(j)))];
        pvh = [pvh, mean(yhatm(q&xdm.time == times(j)))];
    end
    
    subplot(1,2,1)
    hold on
    plot(times, pv, '-o', 'Color', geno2color(genos{i}), 'LineWidth', 3, 'MarkerFaceColor', geno2color(genos{i}), 'MarkerSize', 10);
    
    subplot(1,2,2);
    hold on
    plot(times, pvh, '-o', 'Color', geno2color(genos{i}), 'LineWidth', 3, 'MarkerFaceColor', geno2color(genos{i}), 'MarkerSize', 10);
end

sf = get_standard_figure_font_sizes;
subplot(1,2,1)
axis square
axis tight;
buffer_axis
set(gca, 'XTick', [0 1 5 8]);
set(gca, 'FontSize', sf.axis_tick_labels);
xlabel('Time', 'FontSize', sf.axis_labels);
ylabel(ylab, 'FontSize', sf.axis_labels);
title('Actual', 'FontSize', sf.axis_labels + 2);
v1 = axis;

subplot(1,2,2)
axis square
axis tight;
buffer_axis
set(gca, 'XTick', [0 1 5 8]);
set(gca, 'FontSize', sf.axis_tick_labels);
set(gca, 'YTick', []);
xlabel('Time', 'FontSize', sf.axis_labels);

title('Predicted', 'FontSize', sf.axis_labels + 2);
v2 = axis;

% Adjust axes
subplot(1,2,1)
axis([0 8 min([v1(3) v2(3)]) max([v1(4) v2(4)])]) ;
sub_pos = get(gca,'position'); % get subplot axis position
set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height

subplot(1,2,2)
axis([0 8 min([v1(3) v2(3)]) max([v1(4) v2(4)])]) ;
sub_pos = get(gca,'position'); % get subplot axis position
set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height












end


function f = pred_v_actual_density_plot(ytrue, yhat, varargin)
    % Standardize first so that all genes can be plotted in the same
    % scale.
    tdensity = setParam(varargin, 'transform_density', @(x) x);
    [~, hn] = unix('hostname');
    if strcmpi(strtrim(hn), 'ygritte')
        fontadd = 14;
    else
        fontadd = 0;
    end
    
    
    subsample = setParam(varargin, 'subsample_samples', []);
    if ~isempty(subsample)
        rng('default');
        subi = randsample(size(ytrue,1), subsample, false);
        ytrue = ytrue(subi,:);
        yhat = yhat(subi,:);
    end
    
%     ycat = standardize([ytrue; yhat]);
% 
%     ytrues = ycat(1:size(ytrue,1),:);
%     ycat(1:size(ytrue,1),:) = [];
%     yhats = ycat;

    ytrues = standardize(ytrue);
    yhats = standardize(yhat);


    ytrue_lin = ytrues(:); %reshape(ytrues, size(ytrues,1)*size(ytrues,2), 1);
    yhat_lin = yhats(:); %reshape(yhats, size(yhats,1)*size(yhats,2), 1);
    avgPCC = corr(ytrue_lin, yhat_lin);
    [~,density,X,Y]=kde2d([yhat_lin, ytrue_lin]);
    %density = tdensity(density); % transformation
    
    f = figure;
    sf = get_standard_figure_font_sizes;
    rdensity = tdensity( (density - min(min(density)))/max(max(density)) );
    [hc, hc] = contourf(X,Y,rdensity,100);
    set(hc, 'LineStyle', 'none')
    %colorbar
    %caxis([0 1]);
    axis([-3 3 -3 3])
    axis square;
    colormap(flipud(hot)); %
    %colormap(cbrewer('seq', 'YlOrRd', 100, 'cubic') )
    set(gca, 'XTick', [-3:3]);
    set(gca, 'YTick', -3:3);
    xlabel('Predicted expression', 'FontSize', sf.axis_labels + fontadd);
    ylabel('Actual expression', 'FontSize', sf.axis_labels + fontadd);
    set(gca, 'FontSize', sf.axis_tick_labels + fontadd);
    box off

    txt = text(-2.8, 2.5, sprintf('Avg. PCC = %0.2f', avgPCC));
    set(txt, 'FontSize', sf.axis_tick_labels+6 + fontadd);
    set(txt, 'Color', 'k');

end
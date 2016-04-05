function f = pred_v_actual_density_plot(ytrue, yhat, varargin)
    % Standardize first so that all genes can be plotted in the same
    % scale.
    
    subsample = setParam(varargin, 'subsample_genes', []);
    if ~isempty(subsample)
        rng('default');
        subi = randsample(size(ytrue,2), subsample);
        ytrue = ytrue(:,subi);
        yhat = yhat(:,subi);
    end
    
%     ycat = standardize([ytrue; yhat]);
% 
%     ytrues = ycat(1:size(ytrue,1),:);
%     ycat(1:size(ytrue,1),:) = [];
%     yhats = ycat;

    ytrues = ytrue;
    yhats = yhat;


    ytrue_lin = reshape(ytrues, size(ytrues,1)*size(ytrues,2), 1);
    yhat_lin = reshape(yhats, size(yhats,1)*size(yhats,2), 1);
    avgPCC = corr(ytrue_lin, yhat_lin);
    [~,density,X,Y]=kde2d([yhat_lin, ytrue_lin]);
    
    f = figure;
    sf = get_standard_figure_font_sizes;
    rdensity = (density - min(min(density)))/max(max(density));
    [hc, hc] = contourf(X,Y,rdensity,100);
    set(hc, 'LineStyle', 'none')
    colorbar
    caxis([0 1]);
    axis([-3 3 -3 3])
    axis square;
    colormap(prgn); % cbrewer('seq', 'Purples', 100, 'cubic') 
    xlabel('Predicted expression', 'FontSize', sf.axis_labels);
    ylabel('Actual expression', 'FontSize', sf.axis_labels);
    set(gca, 'FontSize', sf.axis_tick_labels);

    txt = text(-2.8, 2.5, sprintf('Avg. PCC = %0.2f', avgPCC));
    set(txt, 'FontSize', sf.axis_tick_labels+6);
    set(txt, 'Color', 'w');

end
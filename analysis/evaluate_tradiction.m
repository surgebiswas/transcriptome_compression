function [ perfstats ] = evaluate_tradiction( ytrue, yhat, varargin )

    submission_ids = setParam(varargin, 'submissionids', []);
    predids = setParam(varargin, 'keeppredforidx', []);
    calcpredvactualrsq = setParam(varargin, 'calcpredvactualrsq', true);
    predvactualdensityplot = setParam(varargin, 'predvactualdensityplot', true);
    predvactualhistogram = setParam(varargin, 'predvactualhistogram', true);

    if ~isempty(submission_ids)
        fprintf('Calculating intra-submission adjustment.\n');
        % Evaluate submission adjusted (intra-submission) R^2 and slope
        yhat_sa = subadjust(yhat, submission_ids);
        ytrue_sa = subadjust(ytrue, submission_ids);
    end
    
    
    if calcpredvactualrsq
        % Evaluate global R^2.
        fprintf('Evaluating global PCC.\n');
        perfstats.global_Rsq = rsq_and_slope(ytrue,yhat);
        
        if ~isempty(submission_ids)
            % Evaluate intra-submission R^2
            fprintf('Evaluating intra-submission PCC.\n');
            perfstats.sub_adj_Rsq = rsq_and_slope(ytrue_sa, yhat_sa);
        end
        
        if predvactualhistogram
            perfstats.figs.fh_hist_global = pred_v_actual_PCC_histogram(perfstats.global_Rsq);
            glax = axis;
            
            if ~isempty(submission_ids)
                perfstats.figs.fh_hist_subadj = pred_v_actual_PCC_histogram(perfstats.sub_adj_Rsq);
                saax = axis;
                
                set(0, 'currentfigure', perfstats.figs.fh_hist_global);
                axis([0 1 0 max(glax(4), saax(4))]);
                set(0, 'currentfigure', perfstats.figs.fh_hist_subadj);
                axis([0 1 0 max(glax(4), saax(4))]);
                
            end
        end
        

        if ~isempty(predids)
            perfstats.kept_pred_idx = predids;
            perfstats.kept_yhat = yhat(:,predids);
            perfstats.kept_ytrue = ytrue(:,predids);
            % User can submission adjust these on their own.
        end
    end
    
    if predvactualdensityplot
        
        fprintf('Generating density plot for global predictions.\n');
        perfstats.figs.fh_density_global = pred_v_actual_density_plot(ytrue, yhat);
        title('Global');
        glc = caxis;
        cnew = glc;
        
        if ~isempty(submission_ids)
            fprintf('Generating density plot for intra-submission predictions.\n');
            perfstats.figs.fh_density_subadj = pred_v_actual_density_plot(ytrue_sa, yhat_sa);
            title('Intra-submission');
            sac = caxis;
            
            % Need to match the color axes. 
            cnew = min([glc;sac]);
            
            set(0, 'currentfigure', perfstats.figs.fh_density_global);
            caxis(cnew);
            set(0, 'currentfigure', perfstats.figs.fh_density_subadj);
            caxis(cnew);
            
        end
        
        perfstats.figs.cbar = figure;
        colorbar;
        colormap(prgn);
        axis off;
        caxis(cnew);
        sfs = get_standard_figure_font_sizes;
        set(gca, 'FontSize', sfs.axis_tick_labels);
        
        
    end
    
    function f = pred_v_actual_PCC_histogram(pcc)
        f = figure;
        sf = get_standard_figure_font_sizes;
        [N, B] = hist(pcc, 100);
        hBar = bar(B,N,'hist');  
        axis square;
        v = axis;
        axis([0 1 v(3) v(4)]);
        
        pct = prctile(pcc, [10 50 90]);
        
        ind = false(1,length(B));
        for i = 1 : length(pct)
            ind = ind | abs(B-pct(i)) < diff(B(1:2))/2;  %# Find the index of the containing bin
        end
        
        colors = zeros(length(B),3);
        colors(ind,:) = repmat([1 0 1], sum(ind), 1);
        set(hBar, 'FaceVertexCData', colors);
        
        
        pctrnd = round(pct*100)/100;
        
        set(gca, 'XTick', [0 0.25 pctrnd 1]);
%        rotateXLabels(gca, 45);
        set(gca, 'FontSize', sf.axis_tick_labels);
        xlabel('PCC', 'FontSize', sf.axis_labels);
        ylabel('Num. genes', 'FontSize', sf.axis_labels);
    end
    
    
    function f = pred_v_actual_density_plot(ytrue, yhat)
        % Standardize first so that all genes can be plotted in the same
        % scale.
        ycat = standardize([ytrue; yhat]);
        
        ytrues = ycat(1:size(ytrue,1),:);
        ycat(1:size(ytrue,1),:) = [];
        yhats = ycat;
        
        
        ytrue_lin = reshape(ytrues, size(ytrues,1)*size(ytrues,2), 1);
        yhat_lin = reshape(yhats, size(yhats,1)*size(yhats,2), 1);
        avgPCC = corr(ytrue_lin, yhat_lin);
        [~,density,X,Y]=kde2d([yhat_lin, ytrue_lin]);
        
        f = figure;
        sf = get_standard_figure_font_sizes;
        [hc, hc] = contourf(X,Y,density,100);
        set(hc, 'LineStyle', 'none')
        axis([-3 3 -3 3])
        axis square;
        colormap(prgn); % cbrewer('seq', 'Purples', 100, 'cubic') 
        xlabel('Predicted expression', 'FontSize', sf.axis_labels);
        ylabel('Actual expression', 'FontSize', sf.axis_labels);
        set(gca, 'FontSize', sf.axis_tick_labels);
        
        txt = text(-2.8, 2.5, sprintf('Avg. PCC = %0.2f', avgPCC));
        set(txt, 'FontSize', sf.axis_tick_labels);
        set(txt, 'Color', 'w');
        
    end
    

    function [rsq] = rsq_and_slope(ytrue, yhat)
        rsq = zeros(1, size(ytrue,2));
       % sl = rsq;
        for i = 1 : length(rsq)
            rsq(i) = corr(yhat(:,i), ytrue(:,i));
%             c = pca([ytrue(:,i), yhat(:,i)]);
%             sl(i) = c(2,1)/c(1,1);
        end
    end

end


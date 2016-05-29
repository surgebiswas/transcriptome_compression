function  analyze_prosperf3_results( results, qt )
    
    % Genes first.
    [ta, sn] = subadjust(results.target_gene, qt.Submission);
    [pa, ~] = subadjust(results.pred_gene, qt.Submission);
    
    if size(results.target_gene,1) > 3000
        pred_v_actual_density_plot(ta(sn>=2,:), pa(sn>=2,:), 'subsample_samples', 3000);
    else
        pred_v_actual_density_plot(ta(sn>=2,:), pa(sn>=2,:));
    end
    plotSave('figures/prospective_performance/genes_density.png');
    close
    
    
    % Processes next
    [ta, sn] = subadjust(results.target_proc, qt.Submission);
    [pa, ~] = subadjust(results.pred_proc, qt.Submission);
    
    pred_v_actual_density_plot(ta(sn>=2,:), pa(sn>=2,:));
    plotSave('figures/prospective_performance/process_density.png');
    close
    
    % heatmap of processes
    [ri,ci] = hclust(results.target_proc);
    
    imagesc(standardize(results.target_proc(ri,ci)), [-2 2]); colormap(parula);
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca, 'TickLength', [0 0]);
    daspect([1 30*(size(results.target_proc,1)/2597)*(150/size(results.target_proc,2)) 1]);
    plotSave('figures/prospective_performance/process_heatmap_true.png')
    close
    
     
    imagesc(standardize(results.pred_proc(ri,ci)), [-2 2]); colormap(parula);
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    set(gca, 'TickLength', [0 0]);
    daspect([1 30*(size(results.target_proc,1)/2597)*(150/size(results.target_proc,2)) 1]);
    plotSave('figures/prospective_performance/process_heatmap_pred.png')
    close
    
    
    
    


end


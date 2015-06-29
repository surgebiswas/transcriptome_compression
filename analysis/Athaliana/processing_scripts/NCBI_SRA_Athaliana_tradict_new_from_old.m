function NCBI_SRA_Athaliana_tradict_new_from_old( lY,qt,somp, latest_date )

    LAMBDASCHEDULE = logspace(0,2,20);
    
    date_cutoff = datenum(latest_date);  
    dn = datenum(qt.release_date);
    q = dn < date_cutoff;
    
    
    
    
    if false;
        x = lY(q,somp.S);
        model = tratrain(lY(q,:), x, 'lambda', LAMBDASCHEDULE);
        save('NCBI_SRA_Athaliana_tradict_new_from_old.mat', 'model');
    else
        load('NCBI_SRA_Athaliana_tradict_new_from_old.mat');
    end
    
    xval = lY(~q, somp.S);
    yhat = tradict(xval, model);
    ytrue = lY(~q,:);
    
    % Global R^2 over the validation set.
    glrsq = zeros(1, size(ytrue,2));
    for i = 1 : length(glrsq)
        glrsq(i) = corr(yhat(:,i), ytrue(:,i));
    end

    % Submission adjusted R^2
    yhat_sa = subadjust(yhat, qt.Submission(~q));
    ytrue_sa = subadjust(ytrue, qt.Submission(~q));
    sarsq = zeros(1,size(ytrue,2));
    for i = 1 : length(sarsq)
        sarsq(i) = corr( yhat_sa(:,i), ytrue_sa(:,i) );
    end
    
    figure; 
    hold on
    hist(glrsq,80);
    hist(sarsq,80);
    h = findobj(gca,'Type','patch');
    set(h(1), 'FaceColor', 'c', 'EdgeColor', 'c', 'FaceAlpha', 0.6);
    set(h(2), 'FaceColor', 'm', 'EdgeColor', 'm', 'FaceAlpha', 0.6);
    
    box on
    axis square;
    l = legend('Global', 'Intra-submission', 'Location', 'NorthWest');
    set(l, 'Box', 'off');
    set(gca, 'FontSize', 14);
    ylabel('Number of genes', 'FontSize', 16);
    xlabel('R^2', 'FontSize', 16);
    plotSave('figures/tradiction/tradicting_new_from_old_global_and_submission_adjusted_Rsq_histogram.png');
 
   


end


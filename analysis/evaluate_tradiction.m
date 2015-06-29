function [ perfstats ] = evaluate_tradiction( ytrue, yhat, varargin )

    submission_ids = setParam(varargin, 'submissionids', []);
    predids = setParam(varargin, 'keeppredforidx', []);

    % Evaluate global R^2 and global slope.
    [glrsq] = rsq_and_slope(ytrue,yhat);
    perfstats.global_Rsq = glrsq;
    %perfstats.global_slope = glslope;
    
    if ~isempty(submission_ids)
        % Evaluate submission adjusted (intra-submission) R^2 and slope
        yhat_sa = subadjust(yhat, submission_ids);
        ytrue_sa = subadjust(ytrue, submission_ids);
        [sarsq] = rsq_and_slope(ytrue_sa, yhat_sa);
        perfstats.sub_adj_Rsq = sarsq;
        %perfstats.sub_adj_slope = saslope;
    end
    
    
    if ~isempty(predids)
        perfstats.kept_pred_idx = predids;
        perfstats.kept_yhat = yhat(:,predids);
        perfstats.kept_ytrue = ytrue(:,predids);
        % User can submission adjust these on their own.
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


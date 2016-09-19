function [rsq, sl] = rsq_and_slope(ytrue, yhat, varargin)
    
    ab = setParam(varargin, 'abs', false);

    rsq = zeros(1, size(ytrue,2));
    sl = rsq;
    for i = 1 : length(rsq)
        if ab
            rsq(i) = abs(corr(yhat(:,i), ytrue(:,i)));
        else
            rsq(i) = corr(yhat(:,i), ytrue(:,i));
        end

        s = std(ytrue(:,i));
        mask = abs(ytrue(:,i))>s;

        sl(i) = mean(ytrue(mask,i)./yhat(mask,i));
%             c = pca([ytrue(:,i), yhat(:,i)]);
%             sl(i) = c(2,1)/c(1,1);
    end
end
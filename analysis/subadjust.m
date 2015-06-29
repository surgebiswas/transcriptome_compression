 function sy = subadjust(y, sub)
        subu = unique(sub);
        sy = zeros(size(y));
        for j = 1 : length(subu)
            kk = strcmpi(sub, subu(j));
            sy(kk,:) = y(kk,:) - repmat(mean(y(kk,:),1), sum(kk), 1);
        end
    end

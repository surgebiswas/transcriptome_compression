function [ x, xi ] = hiersort( x )
% hierarchically sorts X by first sorting elements in the first column
% then, by sorting the next column where the previous is constant etc.
% works only on cell arrays at the moment.

[x, xi] = hsort(x);


    function [m, mi] = hsort(m)
        [~, mi] = sort(m(:,1));
        m = m(mi,:);
        
        if size(m,2) > 1
            subblock = sort(unique(m(:,1)));
            for i = 1 : length(subblock)
                sbi = strcmpi(m(:,1), subblock{i});
                [a, si] = hsort(m(sbi,2:end));
                m(sbi,:) = [m(sbi,1), a];
                ta = mi(sbi);
                mi(sbi) = ta(si);
            end
        end
    end


end


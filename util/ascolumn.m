function [ x ] = ascolumn( x )
if ~iscolumn(x)
    x = x';
end

end


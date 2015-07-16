function buffer_axis(varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

buffamt = setParam(varargin, 'buffer', 0.05);

v = axis;
if length(v) == 4
    va = abs(v)*buffamt.*[-1 1 -1 1];
elseif length(v) == 6
    va = abs(v)*buffamt.*[-1 1 -1 1 -1 1];
end

axis(v + va);







end


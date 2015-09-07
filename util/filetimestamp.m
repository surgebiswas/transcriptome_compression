function ts  = filetimestamp
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
ts = strrep(strrep(strrep(datestr(now), ' ', '_'), ':', ''), '-', '');

end


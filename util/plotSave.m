function plotSave(outfile, varargin)

% change = ~isstrprop(outfile, 'alphanum');
% outfile(change) = '_';

dpi = setParam(varargin, 'dpi', 300);
g = setParam(varargin, 'fighandle', gcf);

set(g,'units','normalized','outerposition',[0 0 1 1])
print(g, outfile, '-dpng', ['-r', num2str(dpi)]);

end


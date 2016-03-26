function [h, key] = groupPlot( x, y, d, varargin)
% Plots points in x and y according to variable values in design matrix dataset d.

assert(size(x,1) == size(y,1));
assert(size(y,1) == size(d,1));
n = size(x,1);

shapes = {'o', 's', 'v', '+'};
colors = setParam(varargin, 'usecolors', {'r', 'm', 'g', 'b', 'k', 'c', [0.7 0.7 0], [0 0.7 0.7], [0.7 0 0.7], [0.7 0.7 0.7]});
openness = setParam(varargin, 'openvals', {'closed', 'open', 'open_cross'});

linewidth = setParam(varargin, 'linewidth', 2);
msize = setParam(varargin, 'markersize', 10);
shapevar = setParam(varargin, 'shape', []);
colorvar = setParam(varargin, 'color', []);
openvar = setParam(varargin, 'openness', []); % Only for binary categorical variables.
cm = setParam(varargin, 'usecolormap', []); % If passed, flags that colorvar should be 
                                                % treated as numeric and points
                                                % should be colored with its
                                                % associated heat defined in
                                                % 'cm'.


z = setParam(varargin, 'z', []);


if ~isempty(shapevar)
    [shapegroup, key.shapevals]  = categ2design( d(:,{shapevar}), false, 'rmli', false);
else
    shapegroup = [];
end

if ~isempty(cm)
    colors = cm;
    if isnumeric(colors) && size(colors,2) == 3
        % If this is a standard colormap in matrix form, then convert to a cell
        % array with one row, where each column is a separate color of the colormap.
        q = cell(1, size(colors,1));
        for i = 1 : size(colors,1)
            q{i} = colors(i,:);
        end
        colors = q;
    end
    
    cv = eval(['d.', colorvar]);
    delta = (1e-8)*max(cv);
    assert(isnumeric(cv), 'If colormap is passed, then the ''color'' variable must be numeric.');
    
    colorgroup = ceil(size(cm,1)*( (cv + delta) - min(cv))/(max(cv) + delta - min(cv))); % indexes into the colormap;
    key.colorvals = [min(cv), max(cv)];
else
    if ~isempty(colorvar)
        [colorgroup, key.colorvals] = categ2design( d(:,{colorvar}), false, 'rmli', false);
    else
        colorgroup = [];
    end
end

if ~isempty(openvar)
    [opengroup, key.opennessvals] = categ2design( d(:,{openvar}), false, 'rmli', false);
else
    opengroup = [];
end

if isempty(shapegroup)
    shapegroup = true(n,1);
end

if isempty(colorgroup)
    colorgroup = true(n,1);
end

if isempty(opengroup)
    opengroup = true(n,1);
end

assert(size(shapegroup,2) <= length(shapes));
if isempty(cm); assert(size(colorgroup,2) <= length(colors)); end
%assert(size(opengroup,2) <= 2); % can only be open or closed

key.shapes = shapes(1:size(shapegroup,2));
key.openness = openness(1:size(opengroup,2));
if isempty(cm)
    key.colors = colors(1:size(colorgroup,2));
else
    key.colors = cm;
end


h = gca;
hold on
if isempty(cm)
    
    
    for i = 1 : size(shapegroup,2)
        for j = 1 : size(colorgroup,2)
            for k = 1 : size(opengroup,2);
                s = shapegroup(:,i) & colorgroup(:,j) & opengroup(:,k); 
                shape = shapes{i};
                color = colors{j};

                if k == 1 % plot closed (by default)
                    if isempty(z)
                        plot(x(s), y(s), shape, 'Color', color, 'MarkerSize', msize, 'LineWidth', linewidth, 'MarkerFaceColor', color)
                    else
                        plot3(x(s), y(s), z(s), shape, 'Color', color, 'MarkerSize', msize, 'LineWidth', linewidth, 'MarkerFaceColor', color)
                    end

                elseif k == 2 % plot open
                    if isempty(z)
                        plot(x(s), y(s), shape, 'Color', color, 'LineWidth', linewidth, 'MarkerSize', msize)
                    else
                        plot3(x(s), y(s), z(s), shape, 'Color', color, 'LineWidth', linewidth, 'MarkerSize', msize)
                    end
                elseif k == 3 % plot open_cross
                    if isempty(z)
                        plot(x(s), y(s), shape, 'Color', color, 'LineWidth', linewidth, 'MarkerSize', msize)
                        plot(x(s), y(s), '+', 'Color', color, 'LineWidth', linewidth, 'MarkerSize', msize)
                    else
                        plot3(x(s), y(s), z(s), shape, 'Color', color, 'LineWidth', linewidth, 'MarkerSize', msize)
                        plot3(x(s), y(s), z(s), '+', 'Color', color, 'LineWidth', linewidth, 'MarkerSize', msize)
                    end
                    
                end

            end
        end
    end
    
    
else
    
    % 3D plotting currently unsupported when cm passed
    shapegroup = logical(shapegroup);
    opengroup = logical(opengroup);
            
            for k = 1 : length(colorgroup) % colorgroup is an index into the colormap
                shape = shapes{shapegroup(k,:)};
                closeval = find(opengroup(k,:));
                color = colors{colorgroup(k)};
                
                openSetting = openness{closeval};
                
                if strcmpi(openSetting, 'closed')
                    plot(x(k), y(k), shape, 'Color', color, 'LineWidth', linewidth, 'MarkerSize', msize, 'MarkerFaceColor', color);
                elseif strcmpi(openSetting, 'open')   % plot open
                    plot(x(k), y(k), shape, 'Color', color, 'LineWidth', linewidth, 'MarkerSize', msize);
                elseif strcmpi(openSetting, 'open_cross')
                    plot(x(k), y(k), shape, 'Color', color, 'LineWidth', linewidth, 'MarkerSize', msize);
                    plot(x(k), y(k), '+', 'Color', color, 'LineWidth', linewidth, 'MarkerSize', msize);
                end
                
            end

    
    
end




end


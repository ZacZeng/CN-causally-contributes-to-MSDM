function h = errorshade(y1, y2, varargin)
%% errorshade(...) fills the specified area in a shaded color
%
% errorshade(x, y1, y2, ...) plots the shaded area y1 / y2 over x
%
% errorshade(y1, y2) plots the shaded area y1 / y2 over 1..n
%
% Additional arguments are:
%
% errorshade(..., c) sets the color of the shaded area. c needs to be a color
% specification that is accepted by colorget(). By default, [0.5 0.5 0.5] is
% used.
%
% errorshade(..., 'ShadeColor', c) specifies the color with which the given
% color is mixed to get the shading effect. Defaults to [1 1 1].
%
% errorshade(..., 'ShadeStrength', s) specifies the strength of the mixing, that
% is, to which degree the shade color is mixed with the given color. Defaults to
% 0.7.
%
% errorshade(..., 'FixedColor', c) specifies the true color of the shaded area
% without mixing. All other color options are ignored if this option is given.
%
% errorshade(..., 'EdgeColor', c) specifies the edge color of the shaded area.
% By default, the edge color equals the fill color.
%
% errorshade(..., 'ShowInLegend', true) prevents hiding the area when generating
% the legend. By default it is set to false.
%
% h = errorshade(...) returns the handle to the shaded area.


%% process arguments
if ~isvector(y1), error('First argument needs to be a vector'); end
n = length(y1);
p = inputParser;
p.FunctionName = 'errorshade';
p.addRequired('y1', @(x) isvector(x) && length(x) == n);
p.addRequired('y2', @(x) isvector(x) && length(x) == n);
p.addOptional('y3', [], @(x) isvector(x) && length(x) == n);
p.addOptional('color', [0.5 0.5 0.5], @(c) colorget(c, 'validate', true));
p.addParamValue('ShadeColor', [1 1 1], @(c) colorget(c, 'validate', true));
p.addParamValue('ShadeStrength', 0.7, @(x) isscalar(x) && isreal(x));
p.addParamValue('FixedColor', [], @(c) colorget(c, 'validate', true));
p.addParamValue('EdgeColor', [], @(c) colorget(c, 'validate', true));
p.addParamValue('ShowInLegend', false, @(x) isscalar(x) && islogical(x));
p.parse(y1, y2, varargin{:});

% differentiate between (x, y1, y2, ...) and (y1, y2, ...)
if isempty(p.Results.y3)
    x = [];
    y1 = p.Results.y1(:);
    y2 = p.Results.y2(:);
else
    x = p.Results.y1(:);
    y1 = p.Results.y2(:);
    y2 = p.Results.y3(:);
end
area_color = colorget(p.Results.color);
shade_color = colorget(p.Results.ShadeColor);
shade_strength = max(0, min(1, p.Results.ShadeStrength));
color = p.Results.FixedColor;
edge_color = p.Results.EdgeColor;
show_in_legend = p.Results.ShowInLegend;

%% color for shaded area and edge
if isempty(color)
    color = (1 - shade_strength) * area_color + shade_strength * shade_color;
else
    color = colorget(color);
end
if isempty(edge_color), edge_color = color;
else edge_color = colorget(edge_color); end

%% ensure that x is in ascending order, and y1 < y2
% first check x
if isempty(x)
    x = (1:n)';
else
    [x t] = sort(x);
    y1 = y1(t);
    y2 = y2(t);
end
% remove eventual NaN/Inf in y1, y2 and x
inp_fin = isfinite(y1) & isfinite(y2) & isfinite(x);
y1 = y1(inp_fin);  y2 = y2(inp_fin);  x = x(inp_fin);
% make sure that y1 < y2
y_flip = y1 >= y2;
y_tmp = y1(y_flip);
y1(y_flip) = y2(y_flip);
y2(y_flip) = y_tmp;

%% plot shaded area
h = fill([x; flipud(x)], [y1; flipud(y2)], color, 'EdgeColor', edge_color);
% avoid it being displayed in the legend
if ~show_in_legend
    hAnnotation = get(h,'Annotation');
    hLegendEntry = get(hAnnotation','LegendInformation');
    set(hLegendEntry,'IconDisplayStyle','off')
end

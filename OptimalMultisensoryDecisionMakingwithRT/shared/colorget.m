function c = colorget(c, varargin)
%% c = colorget(c, ...)
%
% returns the rgb color code for the given color c.
%
% c can be given as
% 'x' - x is a MATLAB standard color {y, m, c, r, g, b, w, k}
% x - x is a color-identifying integer, starting from 0
% [r g b] - the color as rgb, all within 0 and 1
%
% Additional parameters are
%
% colorget(..., 'validate', true) returns true if the given color is a valid
% color and false otherwise.
%
% colorget(..., 'brighten', bperc) mixes the given color with bperc % of white.
%
% colorget(..., 'darken', dperc) mixes the given color with dperc % of black.
%
% colorget(..., 'grey', gperc) mixes the given color with gperc % of grey.

% matlab pre-defined colors
matlab_map = struct('y', [1 1 0], 'm', [1 0 1], 'c', [0 1 1], 'r', [1 0 0], ...
    'g', [0 1 0], 'b', [0 0 1], 'w', [1 1 1], 'k', [0 0 0]);
% home-made easy distinguishable colors
colors = [1 0 0; 1 0.63 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1; 1 0 1; ...
    1 0.45 0.45; 0.7 1 0.22; 0.32 1 0.68; 0.35 0.63 1; 0.63 0.4 1; ...
    1 0.58 0.75];
col_num = length(colors);
% grey intensity
grey_int = 0.5;


%% process arguments
p = inputParser;
p.FunctionName = 'colorget';
p.addRequired('c', @(x) true);
p.addParamValue('validate', false, @(x) isscalar(x) && islogical(x));
p.addParamValue('brighten', 0, @(x) isscalar(x) && isreal(x));
p.addParamValue('darken', 0, @(x) isscalar(x) && isreal(x));
p.addParamValue('grey', 0, @(x) isscalar(x) && isreal(x));
p.parse(c, varargin{:});

c = p.Results.c;
validate = p.Results.validate;
bperc = max(0, min(1, p.Results.brighten));
dperc = max(0, min(1, p.Results.darken));
gperc = max(0, min(1, p.Results.grey));

%% turn c into rgb format
if ischar(c) && length(c) == 1 && isfield(matlab_map, c)
    c = matlab_map.(c);
elseif ~ischar(c) && isreal(c) && isscalar(c)
    c = mod(round(c) - 1, col_num * 2) + 1;
    if c <= col_num, c = colors(c, :);
    else c = 0.5 * colors(c - col_num, :) + 0.5 * grey_int * ones(1, 3);
    end
elseif ~ischar(c) && isreal(c) && isvector(c) && length(c) == 3
    c = max(0, min(1, c(:)'));
else
    % did not recognise format of c
    if validate
        c = false;
        return
    end
    error('Color format not recognised');
end
if validate
    c = true;
    return
end

%% apply modifiers successively
c = (1 - bperc) * c + bperc * ones(1, 3);
c = (1 - dperc) * c;
c = (1 - gperc) * c + gperc * grey_int * ones(1, 3);
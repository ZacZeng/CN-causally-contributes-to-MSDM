function y = bounded_f(f, x, x_min, x_max, bound_type)
%% applies soft boundaries to function f
%
% f needs to be a function handle that can be evaluated at f(x). If the
% given x is outside [x_min, x_max], f(x) will be evaluated at the closest
% bounday, with an extrapolated increase depending on the distance of the
% boundary and bound_type.
%
% The following bound_types are supported:
% - 'linear' (default): linear increase outside [x_min, x_max]
% - 'exp': neg. exponential increase outside [x_min, x_max]
%
% 2013 Feb  1: Jan Drugowitsch, created version 0.1

% distance & violation of to boundaries
min_diff = x_min - x;
max_diff = x - x_max;
i_min = min_diff > 0;
i_max = max_diff > 0;

% apply boundary
if any(i_min) || any(i_max)
    % evaluate function at closest point at boundary
    x_b = x;
    x_b(i_min) = x_min(i_min);
    x_b(i_max) = x_max(i_max);
    y = f(x_b);
    
    % extrapolate using boundary function
    if nargin < 5, bound_type = 'linear'; end
    switch bound_type
        case 'linear'
            y = y + sum(min_diff(i_min)) + sum(max_diff(i_max));
        case 'exp'
            if y < inf
                % avoid inf in exp
                y = min(realmax, y + sum(exp(min_diff(i_min))-1) + ...
                                     sum(exp(max_diff(i_max))-1));
            end
        otherwise
            error('Unsupported boundary type "%s"', bound_type);
    end
    
else
    % not violating boundaries
    y = f(x);
end
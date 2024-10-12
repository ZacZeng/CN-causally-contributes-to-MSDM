function slsample_state = ...
    slsample_init(logpdf, z0, z_min, z_max, width, adjust_width, max_iter)
%% slsample_state = 
%     slsample_init(logpdf, z0, z_min, z_max, width, adjust_width, max_iter)
%
% initialises the state of the slice sampler
%
% z0 is the starting state (potentially a vector). z0_min and z0_max specify the
% lower and upper bound for z and default to -Inf and +Inf. width is the initial
% slice width and defaults to 10. If adjust_width is true, then the width is
% adjusted between the different sampling steps. max_iter is the maximum number
% of step-out and shinkage steps when adjusting the slice width and defaults to
% 500.

% process arguments and set defaults
if nargin < 7
    max_iter = 500;
    if nargin < 6
        adjust_width = false;
        if nargin < 5
            width = 10;
            if nargin < 4
                z_min = -Inf;
                if nargin < 3
                    z_max = Inf;
                    if nargin < 2
                        error('At least two arguments required');
                    end
                end
            end
        end
    end
end

% turn everytrhing into a column vector
z0 = z0(:);  z_min = z_min(:);  z_max = z_max(:);
width = width(:);

% check z0 and logpdf function
if ~isvector(z0)
    error('z0 needs to be a scalar or vector');
end
if ~isa(logpdf, 'function_handle')
    error('logpdf needs to be a function handle');
end

% replicate elements if z is multidimensional
if isscalar(z_min), z_min = repmat(z_min, size(z0)); end
if isscalar(z_max), z_max = repmat(z_max, size(z0)); end
if isscalar(width), width = repmat(width, size(z0)); end

% generate state structure
slsample_state = struct('logpdf', logpdf, 'logpdf_z', NaN, ...
    'z', z0, 'z_min', z_min, 'z_max', z_max, ...
    'width', width, 'adjust_width', adjust_width, 'max_iter', max_iter);

% add width buffer if width is to be adjusted
if adjust_width
    slsample_state.width_buffer = repmat(width(:), 1, 10);
    if length(z0) > 1
        error('Cannot adjust witdth for multivariate distributions');
    end
end
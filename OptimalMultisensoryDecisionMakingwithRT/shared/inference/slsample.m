function [z_samples sls nevals logpdf_data] = slsample(sls, nsamples, params)
%% [z_samples sls nevals logpdf_data] = slsample(slsample_state, nsamples, params)
%
% returns nsamples (or 1) samples of z, using slice sampling.
%
% sls is the state of the slice sampler and needs to be initialised using
% slsample_init. nsamples is the number of samples drawn, and default to 1 if
% not given. Irrsepective of the shape of z, z_samples is always a matrix with
% the samples given as rows. nevals contains the number of calls to logpdf. If
% params is given, it is given as the second argument to the logpdf function.
%
% sls is the updated state of the sampler, and nevals is the total number of
% calls to the logpdf function. logpdf(z) for the returned sample is accessible
% by sls.logpdf_z.
%
% If logpdf_data is requested, logpdf needs to return two output variables, and
% logpdf_data returns the second output variable from the last call of logpdf.
% This might be useful to facilitate data that has been generated while
% evaluating the log-probability.
%
% slice sampling is performed as described in Neal (2003), using step-out for
% unimodal pdf's and standard shrinkage in all cases.

% process arguments
if nargin < 2
    nsamples = 1;
    if nargin < 1
        error('Too few arguments');
    end
end
return_data = nargout >= 4;

% copy into local namespace for more convenient access
logpdf = sls.logpdf;
logpdf_z = sls.logpdf_z;
z = sls.z;
z_min = sls.z_min;
z_max = sls.z_max;
max_iter = sls.max_iter;
width = sls.width;
adjust_width = sls.adjust_width;
if adjust_width, width_buffer = sls.width_buffer; end

% include additional parameters, if given
if nargin == 3
    logpdf = @(z) logpdf(z, params);
    % additional parameter invalidates previous function evaluation
    logpdf_z = NaN;
end

nevals = 0;
dims = length(z);

% find logpdf at first z, if not yet evaluated
if isnan(logpdf_z)
    logpdf_z = logpdf(z);
    nevals = nevals + 1;
end

% iterate over all required samples
z_samples = zeros(nsamples, dims);
for i = 1:nsamples
    % draw sample exp(u) uniformly from (0,pdf(z))
    u = logpdf_z + log(rand(1, 1));

    % get random interval around z
    z_lower = z - rand(dims, 1) .* width;
    z_upper = z_lower + width;

    % apply limits
    z_lower = max(z_min, z_lower);
    z_upper = min(z_max, z_upper);

    % perform step-out procedure for unimodal pdf's
    if dims == 1

        % only step out lower bound if above lower limit
        if z_lower > z_min
            iter = 0;
            pzl = logpdf(z_lower);
            % continue while logpdf(z_lower) > u
            while iter <= max_iter && pzl > u && z_lower > z_min;
                z_lower = max(z_min, z_lower - width);
                if z_lower > z_min
                    iter = iter + 1;
                    pzl = logpdf(z_lower);
                end
            end
            if iter > max_iter
                error('Step-out of lower bound failed');
            end
            nevals = nevals + iter + 1;
        end

        % only step out upper bound if below upper limit
        if z_upper < z_max
            iter = 0;
            pzu = logpdf(z_upper);
            % continue while logpdf(z_upper) > u
            while iter <= max_iter && pzu > u && z_upper < z_max
                z_upper = min(z_max, z_upper + width);
                if z_upper < z_max
                    iter = iter + 1;
                    pzu = logpdf(z_upper);
                end
            end
            if iter > max_iter
                error('Step-out of upper bound failed');
            end
            nevals = nevals + iter + 1;
        end
    end

    % sample uniformly from [z_lower, z_upper) to get proposal z
    zp = z_lower + rand(dims, 1) .* (z_upper - z_lower);
    if return_data, [logpdf_zp logpdf_data] = logpdf(zp);
    else logpdf_zp = logpdf(zp); end

    % perform shrinkage as long as logpdf(proporal z) < u
    iter = 0;
    while iter <= max_iter && logpdf_zp < u
        % shrink to proposal z
        upper_idx = zp > z;
        z_upper(upper_idx) = zp(upper_idx);
        lower_idx = ~upper_idx;
        z_lower(lower_idx) = zp(lower_idx);
        % sample new proposal z and evaluate logpdf(proposal z)
        zp = z_lower + rand(dims, 1) .* (z_upper - z_lower);
        if return_data, [logpdf_zp logpdf_data] = logpdf(zp);
        else logpdf_zp = logpdf(zp); end
        iter = iter + 1;
    end
    if iter > max_iter
        error('Shrink-in failed');
    end
    nevals = nevals + iter + 1;

    % update z and store
    z = zp;
    logpdf_z = logpdf_zp;
    z_samples(i, :) = z;
    
    % update witdh buffer and width, if used
    if adjust_width
        width_buffer(:,1:end-1) = width_buffer(:,2:end);
        width_buffer(:,end) = z_upper- z_lower;
        width = mean(width_buffer, 2);
    end
end

% update state
sls.z = z;
sls.logpdf_z = logpdf_z;
if adjust_width
    sls.width = width;
    sls.width_buffer = width_buffer;
end
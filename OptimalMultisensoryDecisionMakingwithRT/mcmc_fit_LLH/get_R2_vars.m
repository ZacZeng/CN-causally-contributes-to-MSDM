function [rt_var, pr_var] = get_R2_vars(dstats)
%% returns the variance of psychometric and chronometric data
%
% dstats are the data statistics as returned by get_data_stats

%% compute means
rt_mean = 0;
pr_mean = 0;
if ~isempty(dstats.vis)
    for c_idx = 1:length(dstats.vis)
        rt_mean = rt_mean + ...
            sum(dstats.vis(c_idx).n_corr.*dstats.vis(c_idx).rt_corr) + ...
            sum(dstats.vis(c_idx).n_incorr.*dstats.vis(c_idx).rt_incorr);
        pr_mean = pr_mean + sum(dstats.vis(c_idx).n.*dstats.vis(c_idx).pr);
    end
end
if ~isempty(dstats.vest)
    rt_mean = rt_mean + ...
        sum(dstats.vest.n_corr.*dstats.vest.rt_corr) + ...
        sum(dstats.vest.n_incorr.*dstats.vest.rt_incorr);
    pr_mean = pr_mean + sum(dstats.vest.n.*dstats.vest.pr);
end
if ~isempty(dstats.comb)
    for c_idx = 1:length(dstats.comb)
        rt_mean = rt_mean + ...
            sum(dstats.comb(c_idx).n_corr.*dstats.comb(c_idx).rt_corr) + ...
            sum(dstats.comb(c_idx).n_incorr.*dstats.comb(c_idx).rt_incorr);
        pr_mean = pr_mean + sum(dstats.comb(c_idx).n.*dstats.comb(c_idx).pr);
    end
end
n = dstats.n;
rt_mean = rt_mean / n;
pr_mean = pr_mean / n;


%% compute variances
rt_var = 0;
pr_var = 0;
if ~isempty(dstats.vis)
    for c_idx = 1:length(dstats.vis)
        rt_var = rt_var + ...
            sum(dstats.vis(c_idx).n_corr.*(dstats.vis(c_idx).rt_corr-rt_mean).^2) + ...
            sum(dstats.vis(c_idx).n_incorr.*(dstats.vis(c_idx).rt_incorr-rt_mean).^2);
        pr_var = pr_var + sum(dstats.vis(c_idx).n.*(dstats.vis(c_idx).pr-pr_mean).^2);
    end
end
if ~isempty(dstats.vest)
    rt_var = rt_var + ...
        sum(dstats.vest.n_corr.*(dstats.vest.rt_corr-rt_mean).^2) + ...
        sum(dstats.vest.n_incorr.*(dstats.vest.rt_incorr-rt_mean).^2);
    pr_var = pr_var + sum(dstats.vest.n.*(dstats.vest.pr-pr_mean).^2);
end
if ~isempty(dstats.comb)
    for c_idx = 1:length(dstats.comb)
        rt_var = rt_var + ...
            sum(dstats.comb(c_idx).n_corr.*(dstats.comb(c_idx).rt_corr-rt_mean).^2) + ...
            sum(dstats.comb(c_idx).n_incorr.*(dstats.comb(c_idx).rt_incorr-rt_mean).^2);
        pr_var = pr_var + sum(dstats.comb(c_idx).n.*(dstats.comb(c_idx).pr-pr_mean).^2);
    end
end
rt_var = rt_var / n;
pr_var = pr_var / n;
function [rt_err, pr_err] = get_R2_errs(ds, ms)
%% returns prediction squared errors for given data and model statistics

% compute squared error
rt_err = 0;
pr_err = 0;
if ~isempty(ms.vis)
    for c_idx = 1:length(ms.vis)
        rt_err = rt_err + ...
            sum(ds.vis(c_idx).n_corr.*(ds.vis(c_idx).rt_corr-ms.vis(c_idx).rt_corr).^2) + ...
            sum(ds.vis(c_idx).n_incorr.*(ds.vis(c_idx).rt_incorr-ms.vis(c_idx).rt_incorr).^2);
        pr_err = pr_err + ...
            sum(ds.vis(c_idx).n.*(ds.vis(c_idx).pr-ms.vis(c_idx).pr).^2);
    end
end
if ~isempty(ms.vest)
    rt_err = rt_err + ...
        sum(ds.vest.n_corr.*(ds.vest.rt_corr-ms.vest.rt_corr).^2) + ...
        sum(ds.vest.n_incorr.*(ds.vest.rt_incorr-ms.vest.rt_incorr).^2);
    pr_err = pr_err + ...
        sum(ds.vest.n.*(ds.vest.pr-ms.vest.pr).^2);
end
if ~isempty(ms.comb)
    for c_idx = 1:length(ms.comb)
        rt_err = rt_err + ...
            sum(ds.comb(c_idx).n_corr.*(ds.comb(c_idx).rt_corr-ms.comb(c_idx).rt_corr).^2) + ...
            sum(ds.comb(c_idx).n_incorr.*(ds.comb(c_idx).rt_incorr-ms.comb(c_idx).rt_incorr).^2);
        pr_err = pr_err + ...
            sum(ds.comb(c_idx).n.*(ds.comb(c_idx).pr-ms.comb(c_idx).pr).^2);
    end
end
rt_err = rt_err / ds.n;
pr_err = pr_err / ds.n;
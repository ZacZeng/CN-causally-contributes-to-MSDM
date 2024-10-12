function l = get_LLH(ds, ms)
%% returns the log-likelihood given data and model structure

l = 0;
if ~isempty(ms.vis)
    for c_idx = 1:length(ms.vis)
        l = l + cond_LLH(ds.vis(c_idx), ms.vis(c_idx));
    end
end
if ~isempty(ms.vest)
    for c_idx = 1:length(ms.vis)
        l = l + cond_LLH(ds.vest, ms.vest);
    end
end
if ~isempty(ms.comb)
    for c_idx = 1:length(ms.comb)
        l = l + cond_LLH(ds.comb(c_idx), ms.comb(c_idx));
    end
end


function l = cond_LLH(ds, ms)
%% returns the log-likelihood for a single condition

% introduced, to avoid the log(l_pr) to cause l=-Inf for very small l_pr
min_llh = -realmax/100;
%min_llh = -Inf;

l_rt_corr = -0.5*(log(2*pi*ms.rt_corr_var./ds.n_corr)+...
    ds.n_corr.*(ms.rt_corr-ds.rt_corr).^2./ms.rt_corr_var);
l_rt_incorr = -0.5*(log(2*pi*ms.rt_incorr_var./ds.n_incorr)+...
    ds.n_incorr.*(ms.rt_incorr-ds.rt_incorr).^2./ms.rt_incorr_var);
% l_rt_corr = -0.5*(log(2*pi*ms.rt_corr_var./ds.n_corr)+...
%     ds.n_corr.*(ms.rt_corr-ds.rt_corr).^2);
% l_rt_incorr = -0.5*(log(2*pi*ms.rt_incorr_var./ds.n_incorr)+...
%     ds.n_incorr.*(ms.rt_incorr-ds.rt_incorr).^2);
l_pr = -0.5*(ds.n .* (ds.pr-ms.pr).^2./ms.rt_corr_var)/2;
% l_pr = log(binopdf(ds.nr,ds.n,ms.pr));
l = sum(max(min_llh,l_rt_corr(ds.n_corr>0)))+...
    sum(max(min_llh,l_rt_incorr(ds.n_incorr>0)))+sum(max(min_llh,l_pr));

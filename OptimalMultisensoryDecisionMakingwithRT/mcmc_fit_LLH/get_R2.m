function r2 = get_R2(subj_id, model_fn)
%% returns R2 measure for given subject / model combination

% load fit data
model_name = func2str(model_fn);
m = load(['fit_data' filesep model_name '_' subj_id '.mat']);
% compute r2
[rt_var, pr_var] = get_R2_vars(m.dstats);
[rt_err, pr_err] = get_R2_errs(m.dstats, m.mstats);
r2 = 1 - 0.5 * (rt_err/rt_var + pr_err/pr_var);

function [ms, biases] = mf_velow(d, p)
%% returns the model predictions for given subject data and parameters
%
% this model impements the predictions that are optimal besides assuming
% that both modalities are weighted by the velocity.
%
% For initialisation, the model assumes to find the fits for
% mf_accw_powcomb in the data folder. It uses those to find adequate
% initialisation parameters per subject.
%
% the model uses the following parameters:
%
% coh_num x k_vis
% coh_num x bound_vis
% k_vest
% bound_vest
% coh_num x bound_comb
% tnd_vis, tnd_vest, tnd_comb
% (2 x coh_num + 1) biases
% pl
%
% resulting in 5 x coh_num + 7 parameters.
%
% d is a structured, expected to contain
% - cohs: sequence of used coherences
% - hs: sequence of headings
% - delta_t: step size used for computations
% - t_max: largest time to compute distributions for
% - v_t: velocity profile in steps of delta_t
% - a_t: acceleration profile in steps of delta_t

% PREV
%Furthermore, we
% assume the k's and bounds to relate to coherence by
% k_vis(c) = a_vis c^gam_vis / sqrt(1+b_vis c^gam_vis)
% bound_vis(c) = g_vis / sqrt(1+b_vis c^gam_vis)
% bound_comb(c) = g_comb / sqrt(1+b_comb c^gam_comb)
%
% the model uses the following parameters:
%
% a_vis (1), b_vis (2), g_vis (3), gam_vis (4)
% k_vest (5), bound_vest (6)
% b_comb (7), g_comb (8), gam_comb (9)
% tnd_vis (10), tnd_vest (11), tnd_comb (12)
% (2 x coh_num + 1) biases
% pl
%
% resulting in 2 x coh_num + 14 parameters.


%% return model parameters if required
if nargin < 2
    cs = d.cohs;
    c_num = length(cs);
    % load best-found parameters from velow_powcomb fit for initialisation
    prev_p = load(['fit_data' filesep 'mf_velow_powcomb_' d.subj_id '.mat']);
    assert(c_num == length(prev_p.cohs));
    prev_p = prev_p.best_p;
    % build initial data structure
    kvis_ini = prev_p(1)*cs.^prev_p(4)./sqrt(1+prev_p(2)*cs.^prev_p(4));
    boundvis_ini = prev_p(3)./sqrt(1+prev_p(2)*cs.^prev_p(4));
    boundcomb_ini = prev_p(8)./sqrt(1+prev_p(7)*cs.^prev_p(9));
    ms.p_min = [(0.00001*ones(1,c_num)) (0.001*ones(1,c_num)) 0.00001 0.001 (0.001*ones(1,c_num))   0   0   0 (-5*ones(1,2*c_num+1))   0];
    ms.p_max = [   (1000*ones(1,c_num))    (10*ones(1,c_num))    1000    10    (10*ones(1,c_num))   1   1   1  (5*ones(1,2*c_num+1)) 0.2];
    ms.p_ini = [kvis_ini boundvis_ini prev_p(5:6) boundcomb_ini prev_p(10:end)];
    ms.p_w   = [          (2+150*cs.^3)  (0.02*ones(1,c_num))      20  0.03  (0.02*ones(1,c_num)) 0.2 0.2 0.2    (ones(1,2*c_num+1)) 0.1];
    ms.p_names = cat(2,...
        arrayfun(@(i) sprintf('kvis%d',i), 1:c_num, 'UniformOutput', false),...
        arrayfun(@(i) sprintf('boundvis%d',i), 1:c_num, 'UniformOutput', false),...
        {'kvest','boundvest'},...
        arrayfun(@(i) sprintf('boundcomb%d',i), 1:c_num, 'UniformOutput', false),...
        {'tndvis','tndvest','tndcomb'},...
        arrayfun(@(i) sprintf('biasvis%d',i), 1:c_num, 'UniformOutput', false),...
        {'biasvest'},...
        arrayfun(@(i) sprintf('biascomb%d',i), 1:c_num, 'UniformOutput', false),...
        {'pl'});
    return
end


%% decode parameters
c = d.cohs;
c_num = length(c);
k_vis = p(1:c_num);
bound_vis = p((c_num+1):(2*c_num));
k_vest = p(2*c_num+1);
bound_vest = p(2*c_num+2);
bound_comb = p((2*c_num+3):(3*c_num+2));
tnds = p((3*c_num+3):(3*c_num+5));
biases = struct('vis', p((3*c_num+6):(4*c_num+5)), 'vest', p(4*c_num+6), ...
    'comb', p((4*c_num+7):(5*c_num+6)));
pl = p(5*c_num+7);


%% compute optimal model predictions
% optimal sensitivity & bounds
k_comb = sqrt(k_vis.^2 + k_vest^2);

% velocity / combined reliability profile
a_vis = d.v_t;
a_vest = d.v_t;
a_comb = zeros(c_num,length(a_vis));
for c_idx = 1:c_num
    a_comb(c_idx,:) = sqrt((k_vis(c_idx)/k_comb(c_idx))^2*a_vis.^2 + ...
        (k_vest/k_comb(c_idx))^2*a_vest.^2);
end

% model predictions per condition / coherence
delta_t = d.delta_t;
t_max = d.t_max;
hs = d.hs;
vis(1,c_num) = struct('pr', [], 'rt_corr', [], 'rt_incorr', [],...
    'rt_corr_var', [], 'rt_incorr_var', []);
comb(1,c_num) = struct('pr', [], 'rt_corr', [], 'rt_incorr', [],...
    'rt_corr_var', [], 'rt_incorr_var', []);
for c_idx = 1:c_num
    vis(c_idx) = model_pred(a_vis,bound_vis(c_idx),delta_t,t_max,tnds(1),pl,...
        k_vis(c_idx)*sin((hs+biases.vis(c_idx))*pi/180));
    comb(c_idx) = model_pred(a_comb(c_idx,:),bound_comb(c_idx),delta_t,t_max,tnds(3),pl,...
        k_comb(c_idx)*sin((hs+biases.comb(c_idx))*pi/180));
end
vest = model_pred(a_vest,bound_vest,delta_t,t_max,tnds(2),pl,...
    k_vest*sin((hs+biases.vest)*pi/180));

ms = struct('vis', vis, 'vest', vest, 'comb', comb);



function m = model_pred(a, bound, delta_t, t_max, tnd, pl, ks)
%% returns model predictions as computed from DM RT distributions

% reserve space
k_num = length(ks);
rt_corr = zeros(1, k_num);    rt_corr_var = zeros(1, k_num);
rt_incorr = zeros(1, k_num);  rt_incorr_var = zeros(1, k_num);
pr = zeros(1, k_num);

% compute per k
for k_idx = 1:k_num
    [g1, g2] = ddm_rt_dist(a, bound, delta_t, t_max, ks(k_idx), 'mnorm', 'yes');
    % compute stats
    g1(isnan(g1) | isinf(g1)) = 0;
    pr(k_idx) = sum(g1)*delta_t;
    g1 = g1 / sum(g1);  g1(isnan(g1) | isinf(g1)) = 0;
    g2 = g2 / sum(g2);  g2(isnan(g2) | isinf(g2)) = 0;
    ts = (1:length(g1))*delta_t - delta_t/2;  ts2 = ts.^2;
    rt_corr(k_idx) = sum(g1.*ts);
    rt_incorr(k_idx) = sum(g2.*ts);
    rt_corr_var(k_idx) = max(1e-20, sum(g1.*ts2) - rt_corr(k_idx)^2);
    rt_incorr_var(k_idx) = max(1e-20, sum(g2.*ts2) - rt_incorr(k_idx)^2);
end

% add lapses and non-decision time
m = struct('pr', (1-pl)*pr + pl*0.5, ...
    'rt_corr', rt_corr+tnd, 'rt_incorr', rt_incorr+tnd, ...
    'rt_corr_var', rt_corr_var, 'rt_incorr_var', rt_incorr_var);
function [ms, biases] = mf_sepw_powsep(d, p)
%% returns the model predictions for given subject data and parameters
%
% this model explicitly fits the weight in the combined condition without
% optimality considerations. I assumes k's and bounds in the unimodal
% condition to relate to coherence by
% k_vis(c) = a_vis c^gam1_vis / sqrt(1+b_vis c^gam2_vis)
% bound_vis(c) = g_vis / sqrt(1+b_vis c^gam2_vis)
% bound_comb(c) = g_comb / sqrt(1+b_comb c^gam_comb)
%
% For initialisation, the model assumes to find the fits for
% mf_sepk_powsep in the data folder. It uses those to find adequate
% initialisation parameters per subject.
%
% the model uses the following parameters:
%
% a_vis (1), b_vis (2), g_vis (3), gam1_vis (4), gam2_vis (5)
% k_vest (6), bound_vest (7)
% w_vis (8 - 7+coh_num)
% b_comb (8+coh_num), g_comb (9+coh_num), gam_comb (10+coh_num)
% tnd_vis, tnd_vest, tnd_comb (11+coh_num - 13+coh_num)
% (2 x coh_num + 1) biases
% pl
%
% resulting in 3 x coh_num + 15 parameters.
%
% d is a structured, expected to contain
% - cohs: sequence of used coherences
% - hs: sequence of headings
% - delta_t: step size used for computations
% - t_max: largest time to compute distributions for
% - v_t: velocity profile in steps of delta_t
% - a_t: acceleration profile in steps of delta_t
% - subj_id: the subject id to fit



%% return model parameters if required
if nargin < 2
    cs = d.cohs(:)';
    c_num = length(cs);
    
    % load best-found parameters from sepk fit for initialisation
    sepk_p = load(['fit_data' filesep 'mf_sepk_powsep_' d.subj_id '.mat']);
    assert(c_num == length(sepk_p.cohs));
    sepk_p = sepk_p.best_p;
    % build initial data structure
    sepk_kvis = sepk_p(1)*cs.^sepk_p(4)./sqrt(1+sepk_p(2)*cs.^sepk_p(5));
    sepk_kvest = sepk_p(6);
    w_vis_ini = sepk_kvis.^2./(sepk_kvis.^2+sepk_kvest^2);
    sepk_p(8:(c_num+7)) = w_vis_ini;
    ms = struct();
    ms.p_min = [0.0001  0 0.001   0   0 0.00001 0.001       zeros(1,c_num)    0 0.001   0   0   0   0 (-5*ones(1,2*c_num+1))   0];
    ms.p_max = [ 10000 50    10  10  10    1000    10        ones(1,c_num)  200    10  10   1   1   1  (5*ones(1,2*c_num+1)) 0.2];
    ms.p_ini = sepk_p;
    ms.p_w   = [    50  1  0.05 0.5 0.5      20  0.03 (0.02*ones(1,c_num))   10  0.05 0.5 0.2 0.2 0.2    (ones(1,2*c_num+1)) 0.1];
    ms.p_names = cat(2,...
        {'avis','bvis','gvis','gam1vis','gam2vis'},...
        {'kvest','boundvest'},...
        arrayfun(@(i) sprintf('wvis%d',i), 1:c_num, 'UniformOutput', false), ...
        {'bcomb','gcomb','gamcomb'},...
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
k_vis = p(1)*c.^p(4)./sqrt(1+p(2)*c.^p(5));
bound_vis = p(3)./sqrt(1+p(2)*c.^p(5));
k_vest = p(6);
bound_vest = p(7);
w_vis = p(8:(7+c_num));
w_vest = 1-w_vis;
bound_comb = p(9+c_num)./sqrt(1+p(8+c_num)*c.^p(10+c_num));
tnds = p((11+c_num):(13+c_num));
biases = struct('vis', p((14+c_num):(2*c_num+13)), 'vest', p(2*c_num+14), ...
    'comb', p((2*c_num+15):(3*c_num+14)));
pl = p(3*c_num+15);


%% compute model predictions
% velocity / acceleration
a_vis = d.v_t(:)';
a_vest = d.a_t(:)';
% noise variance in comb condition, coh_idx x time
sig2_comb = bsxfun(@times,w_vis(:),a_vis.^2) + ...
    bsxfun(@times,w_vest(:),a_vest.^2);


% model predictions per condition / coherence
delta_t = d.delta_t;
t_max = d.t_max;
hs = d.hs(:)';
vis(1,c_num) = struct('pr', [], 'rt_corr', [], 'rt_incorr', [],...
    'rt_corr_var', [], 'rt_incorr_var', []);
comb(1,c_num) = struct('pr', [], 'rt_corr', [], 'rt_incorr', [],...
    'rt_corr_var', [], 'rt_incorr_var', []);
for c_idx = 1:c_num
    vis(c_idx) = model_pred1(a_vis,bound_vis(c_idx),delta_t,t_max,tnds(1),pl,...
        k_vis(c_idx)*sin((hs+biases.vis(c_idx))*pi/180));
    % x_comb = sqrt(w_vis) x_vis + sqrt(w_vest) x_vest, resulting in
    % drift = (sqrt(w_vis) a_vis^2 k_vis + sqrt(w_vest) a_vest^2 k_vest) sin(h+bias)
    % variance = w_vis a_vis^2 + w_vest a_vest^2
    drifts = bsxfun(@times,...
        sqrt(w_vis(c_idx))*a_vis.^2*k_vis(c_idx)+sqrt(w_vest(c_idx))*a_vest.^2*k_vest,...
        sin((hs'+biases.comb(c_idx))*pi/180));
    comb(c_idx) = model_pred2(drifts,sig2_comb(c_idx,:),bound_comb(c_idx),...
        delta_t,t_max,tnds(3),pl);
end
vest = model_pred1(a_vest,bound_vest,delta_t,t_max,tnds(2),pl,...
    k_vest*sin((hs+biases.vest)*pi/180));

ms = struct('vis', vis, 'vest', vest, 'comb', comb);


function m = model_pred1(a, bound, delta_t, t_max, tnd, pl, ks)
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


function m = model_pred2(drifts, sig2, bound, delta_t, t_max, tnd, pl)
%% returns model predictions as computed from DM RT distributions

% reserve space
k_num = size(drifts,1);
rt_corr = zeros(1, k_num);    rt_corr_var = zeros(1, k_num);
rt_incorr = zeros(1, k_num);  rt_incorr_var = zeros(1, k_num);
pr = zeros(1, k_num);

% compute per k
for k_idx = 1:k_num
    [g1, g2] = ddm_rt_dist_full(drifts(k_idx,:),sig2,-bound,bound,0,0,...
        delta_t,t_max,'mnorm','yes');
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
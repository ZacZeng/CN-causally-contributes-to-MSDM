function [ms, biases] = mf_sepk_nobias_onetnd(d, p)
%% returns the model predictions for given subject data and parameters
%
% this model does not assume any relation between unimodal and combined k's
% and bounds. Further, it does not assume any biases, and a single, shared
% non-decision time.
%
% For initialisation, the model assumes to find the fits for mf_sepk in the
% data folder. It uses those to find adequate initialisation parameters per
% subject.
%
%
% the model uses the following parameters:
%
% coh_num x k_vis
% coh_num x bound_vis
% k_vest
% bound_vest
% coh_num x k_comb
% coh_num x bound_comb
% tnd
% pl
%
% resulting in 4 x coh_num + 4 parameters.
%
% d is a structured, expected to contain
% - cohs: sequence of used coherences
% - hs: sequence of headings
% - delta_t: step size used for computations
% - t_max: largest time to compute distributions for
% - v_t: velocity profile in steps of delta_t
% - a_t: acceleration profile in steps of delta_t


%% return model parameters if required
if nargin < 2
    cs = d.cohs;
    c_num = length(cs);
    % load best-found parameters from sepk fit for initialisation
    sepk_p = load(['fit_data' filesep 'mf_sepk_' d.subj_id '.mat']);
    assert(c_num == length(sepk_p.cohs));
    sepk_p = sepk_p.best_p;
    assert(length(sepk_p) == 6*c_num+7);
    ms = struct();
    ms.p_min = [(0.00001*ones(1,c_num)) (0.001*ones(1,c_num)) 0.00001 0.001 (0.00001*ones(1,c_num)) (0.001*ones(1,c_num))   0   0];
    ms.p_max = [   (1000*ones(1,c_num))    (10*ones(1,c_num))    1000    10    (1000*ones(1,c_num))    (10*ones(1,c_num))   1 0.2];
    ms.p_ini = [sepk_p(1:(4*c_num+2)) mean(sepk_p((4*c_num+3):(4*c_num+5))) sepk_p(end)];
    ms.p_w   = [          (2+150*cs.^3)  (0.02*ones(1,c_num))      20  0.03           (2+150*cs.^3)  (0.02*ones(1,c_num)) 0.2 0.1];
    ms.p_names = cat(2,...
        arrayfun(@(i) sprintf('kvis%d',i), 1:c_num, 'UniformOutput', false),...
        arrayfun(@(i) sprintf('boundvis%d',i), 1:c_num, 'UniformOutput', false),...
        {'kvest','boundvest'},...
        arrayfun(@(i) sprintf('kcomb%d',i), 1:c_num, 'UniformOutput', false),...
        arrayfun(@(i) sprintf('boundcomb%d',i), 1:c_num, 'UniformOutput', false),...
        {'tnd','pl'});
    return
end


%% decode parameters
c = d.cohs;
c_num = length(c);
k_vis = p(1:c_num);
bound_vis = p((c_num+1):(2*c_num));
k_vest = p(2*c_num+1);
bound_vest = p(2*c_num+2);
k_comb = p((2*c_num+3):(3*c_num+2));
bound_comb = p((3*c_num+3):(4*c_num+2));
tnd = p(4*c_num+3);
pl = p(4*c_num+4);
biases = struct('vis', zeros(1,c_num), 'vest', 0, ...
    'comb', zeros(1,c_num));



%% compute model predictions
% velocity / acceleration / combined reliability profile
a_vis = d.v_t;
a_vest = d.a_t;
a_comb = zeros(c_num,length(a_vis));
for c_idx = 1:c_num
    a_comb(c_idx,:) = sqrt((k_vis(c_idx)/k_comb(c_idx))^2*a_vis.^2 + ...
        (k_vest/k_comb(c_idx))^2*a_vest.^2);
end

% model predictions per condition / coherence
delta_t = d.delta_t;
t_max = d.t_max;
sin_hs = sin(d.hs*pi/180);
vis(1,c_num) = struct('pr', [], 'rt_corr', [], 'rt_incorr', [],...
    'rt_corr_var', [], 'rt_incorr_var', []);
comb(1,c_num) = struct('pr', [], 'rt_corr', [], 'rt_incorr', [],...
    'rt_corr_var', [], 'rt_incorr_var', []);
for c_idx = 1:c_num
    vis(c_idx) = model_pred(a_vis,bound_vis(c_idx),delta_t,t_max,tnd,pl,...
        k_vis(c_idx)*sin_hs);
    comb(c_idx) = model_pred(a_comb(c_idx,:),bound_comb(c_idx),delta_t,t_max,tnd,pl,...
        k_comb(c_idx)*sin_hs);
end
vest = model_pred(a_vest,bound_vest,delta_t,t_max,tnd,pl,...
    k_vest*sin_hs);

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
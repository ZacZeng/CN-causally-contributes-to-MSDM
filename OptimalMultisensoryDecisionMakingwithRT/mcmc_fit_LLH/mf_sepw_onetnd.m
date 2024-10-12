function [ms, biases] = mf_sepw_onetnd(d, p)
%% returns the model predictions for given subject data and parameters
%
% this model explicitly fits the weight in the combined condition without
% optimality considerations. Furthermore, it models the bounds for each
% condition/coherence separately, and assumes a single, shared non-decision
% time.
%
% For initialisation, the model assumes to find the fits for mf_sepk in the
% data folder. It uses those to find adequate initialisation parameters per
% subject.
%
% the model uses the following parameters:
%
% coh_num x k_vis
% coh_num x bound_vis
% k_vest
% bound_vest
% coh_num x w_vis (weight assigned to visual modality in comb. cond)
% coh_num x bound_comb
% tnd_vis, tnd_vest, tnd_comb
% (2 x coh_num + 1) biases
% pl
%
% resulting in 6 x coh_num + 5 parameters.
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
    cs = d.cohs;
    c_num = length(cs);
    % load best-found parameters from sepk fit for initialisation
    sepk_p = load(['fit_data' filesep 'mf_sepk_' d.subj_id '.mat']);
    assert(c_num == length(sepk_p.cohs));
    sepk_p = sepk_p.best_p;
    % build initial data structure
    w_vis_ini = sepk_p(1:c_num).^2./(sepk_p(1:c_num).^2+sepk_p(2*c_num+1)^2);
    sepk_p((2*c_num+3):(3*c_num+2)) = w_vis_ini;
    ms = struct();
    ms.p_min = [(0.00001*ones(1,c_num)) (0.001*ones(1,c_num)) 0.00001 0.001      zeros(1,c_num) (0.001*ones(1,c_num))   0 (-5*ones(1,2*c_num+1))   0];
    ms.p_max = [   (1000*ones(1,c_num))    (10*ones(1,c_num))    1000    10       ones(1,c_num)    (10*ones(1,c_num))   1  (5*ones(1,2*c_num+1)) 0.2];
    ms.p_ini = [sepk_p(1:(4*c_num+2)) mean(sepk_p((4*c_num+3):(4*c_num+5))) sepk_p((4*c_num+6):end)];
    ms.p_w   = [          (2+150*cs.^3)  (0.02*ones(1,c_num))      20  0.03 (0.2*ones(1,c_num))  (0.02*ones(1,c_num)) 0.2    (ones(1,2*c_num+1)) 0.1];
    ms.p_names = cat(2,...
        arrayfun(@(i) sprintf('kvis%d',i), 1:c_num, 'UniformOutput', false),...
        arrayfun(@(i) sprintf('boundvis%d',i), 1:c_num, 'UniformOutput', false),...
        {'kvest','boundvest'},...
        arrayfun(@(i) sprintf('wvis%d',i), 1:c_num, 'UniformOutput', false),...
        arrayfun(@(i) sprintf('boundcomb%d',i), 1:c_num, 'UniformOutput', false),...
        {'tnd'},...
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
w_vis = p((2*c_num+3):(3*c_num+2));
w_vest = 1-w_vis;
bound_comb = p((3*c_num+3):(4*c_num+2));
tnd = p(4*c_num+3);
biases = struct('vis', p((4*c_num+4):(5*c_num+3)), 'vest', p(5*c_num+4), ...
    'comb', p((5*c_num+5):(6*c_num+4)));
pl = p(6*c_num+5);


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
    vis(c_idx) = model_pred1(a_vis,bound_vis(c_idx),delta_t,t_max,tnd,pl,...
        k_vis(c_idx)*sin((hs+biases.vis(c_idx))*pi/180));
    % x_comb = sqrt(w_vis) x_vis + sqrt(w_vest) x_vest, resulting in
    % drift = (sqrt(w_vis) a_vis^2 k_vis + sqrt(w_vest) a_vest^2 k_vest) sin(h+bias)
    % variance = w_vis a_vis^2 + w_vest a_vest^2
    drifts = bsxfun(@times,...
        sqrt(w_vis(c_idx))*a_vis.^2*k_vis(c_idx)+sqrt(w_vest(c_idx))*a_vest.^2*k_vest,...
        sin((hs'+biases.comb(c_idx))*pi/180));
    comb(c_idx) = model_pred2(drifts,sig2_comb(c_idx,:),bound_comb(c_idx),...
        delta_t,t_max,tnd,pl);
end
vest = model_pred1(a_vest,bound_vest,delta_t,t_max,tnd,pl,...
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
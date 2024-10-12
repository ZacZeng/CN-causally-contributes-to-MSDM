function [ms, biases] = mf_unimodal(d, p)
%% returns the model predictions for given subject data and parameters
%
% this model fits the unimodal conditions only, and fit k's and bounds
% separately for these conditions.
%
% the model uses the following parameters:
%
% coh_num x k_vis
% coh_num x bound_vis
% k_vest
% bound_vest
% tnd_vis, tnd_vest
% (coh_num + 1) biases
% pl
%
% resulting in 3 x coh_num + 6 parameters.
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
    ms = struct();
    ms.p_min = [(0.00001*ones(1,c_num)) (0.001*ones(1,c_num)) 0.00001 0.001   0   0 (-5*ones(1,c_num+1))   0];
    ms.p_max = [   (1000*ones(1,c_num))    (10*ones(1,c_num))    1000    10   1   1  (5*ones(1,c_num+1)) 0.2];
    ms.p_ini = [          (1+500*cs.^2)    (0.25*(1-cs.^2.2))      40  0.25 0.2 0.2     zeros(1,c_num+1)   0];
    ms.p_w   = [          (2+150*cs.^3)  (0.02*ones(1,c_num))      20  0.03 0.2 0.2    (ones(1,c_num+1)) 0.1];
    ms.p_names = cat(2,...
        arrayfun(@(i) sprintf('kvis%d',i), 1:c_num, 'UniformOutput', false),...
        arrayfun(@(i) sprintf('boundvis%d',i), 1:c_num, 'UniformOutput', false),...
        {'kvest','boundvest'},...
        {'tndvis','tndvest'},...
        arrayfun(@(i) sprintf('biasvis%d',i), 1:c_num, 'UniformOutput', false),...
        {'biasvest'},...
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
tnds = p((2*c_num+3):(2*c_num+4));
biases = struct('vis', p((2*c_num+5):(3*c_num+4)), 'vest', p(3*c_num+5), ...
    'comb', []);
pl = p(3*c_num+6);


%% compute optimal model predictions
% velocity / acceleration / combined reliability profile
a_vis = d.v_t;
a_vest = d.a_t;

% model predictions per condition / coherence
delta_t = d.delta_t;
t_max = d.t_max;
hs = d.hs;
vis(1,c_num) = struct('pr', [], 'rt_corr', [], 'rt_incorr', [],...
    'rt_corr_var', [], 'rt_incorr_var', []);
for c_idx = 1:c_num
    vis(c_idx) = model_pred(a_vis,bound_vis(c_idx),delta_t,t_max,tnds(1),pl,...
        k_vis(c_idx)*sin((hs+biases.vis(c_idx))*pi/180));
end
vest = model_pred(a_vest,bound_vest,delta_t,t_max,tnds(2),pl,...
    k_vest*sin((hs+biases.vest)*pi/180));

ms = struct('vis', vis, 'vest', vest, 'comb', []);



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
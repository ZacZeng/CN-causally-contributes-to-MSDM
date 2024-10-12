function dstats = get_data_stats(d, biases)
%% returns the strutured stats that are used for the fits
%
% d is the subject data, as returned by load_condi_data
%
% for each condition, the follow vectors (over heading) are returned
% - hs: headings
% - rt_corr, rt_incorr: mean reaction times, correct/incorrect choices
% - rt_corr_var, rt_incorr_var: corresponding variances
% - pr: probability of rightwards choice
% - n, n_corr, n_incorr: total number of trials, correct/incorrect trials
% - nr: number of trials 'right' option was chosen
%
% If biases are given, they are expected to be a structure with fields
% - vis: array with one element per coherence
% - vest: double
% - comb: array with one element per coherence


%% process arguments
if nargin < 2    
    % biases not given - consider all conditions
    vis_num = length(d.vis);
    comb_num = length(d.comb);
    biases = struct('vest', 0, 'vis', zeros(1, vis_num), ...
        'comb', zeros(1, comb_num));
else
    % biases determine number of conditions
    vis_num = sum(~isnan(biases.vis));
    vest_num = sum(~isnan(biases.vest));
    comb_num = sum(~isnan(biases.comb));
end


%% collect stats per condition / modalitiy
n = 0;
if (vis_num > 0) && ~isempty(biases.vis)
    vis(1, vis_num) = struct('hs', [], ...
        'rt_corr', [], 'rt_incorr', [], 'rt_corr_var', [], 'rt_incorr_var', [], ...
        'pr', [], 'n', [], 'nr', [], 'n_corr', [], 'n_incorr', []);
    for c_idx = 1:vis_num
        vis(c_idx) = compute_cond_stats(d.vis(c_idx), biases.vis(c_idx));
        n = n + sum(vis(c_idx).n);
    end
else
    vis = [];
end
if ~isempty(d.vest) && ~isempty(biases.vest)
    vest = compute_cond_stats(d.vest, biases.vest);
    n = n + sum(vest.n);
else
    vest = [];
end
if (comb_num > 0) && ~isempty(biases.comb)
    comb(1, comb_num) = struct('hs', [], ...
        'rt_corr', [], 'rt_incorr', [], 'rt_corr_var', [], 'rt_incorr_var', [], ...
        'pr', [], 'n', [], 'nr', [], 'n_corr', [], 'n_incorr', []);
    for c_idx = 1:comb_num
        comb(c_idx) = compute_cond_stats(d.comb(c_idx), biases.comb(c_idx));
        n = n + sum(comb(c_idx).n);
    end
else
    comb = [];
end


%% build data structure and return
dstats = struct('vis', vis, 'vest', vest, 'comb', comb, 'n', n);


function cond_stats = compute_cond_stats(d, bias)
%% returns statistics per condition

% split into correct / wrong
es = d.es;
hs = unique(es);
h_num = length(hs);
% % Changed by ZZ @ 20231008 
% % Why this? d.choice is already boolean value indicate correct or not 
% corr_choice = logical(d.choice); 
corr_choice = es >= -bias & logical(d.choice) | es < -bias & ~logical(d.choice);

% sort trials by heading / correct / wrong
rt_corr = zeros(1, h_num);
rt_corr_var = zeros(1, h_num);
rt_incorr = zeros(1, h_num);
rt_incorr_var = zeros(1, h_num);
n = zeros(1, h_num);
nr = zeros(1, h_num);
n_corr = zeros(1, h_num);
for h_idx = 1:h_num
    % choose trials for heading / correct & count them
    h_trials = es == hs(h_idx);
    corr_trials = h_trials & corr_choice;
    n(h_idx) = sum(h_trials);
    n_corr(h_idx) = sum(corr_trials);
    nr(h_idx) = sum(h_trials & d.choice);
    % reaction time, correct trials
    rt = d.rt(corr_trials);
    if ~isempty(rt)
        rt_corr(h_idx) = mean(rt);
        rt_corr_var(h_idx) = var(rt);
    end
    % reaction time, incorrect trials
    rt = d.rt(h_trials & ~corr_choice);
    if ~isempty(rt)
        rt_incorr(h_idx) = mean(rt);
        rt_incorr_var(h_idx) = var(rt);
    end
end
n_incorr = n - n_corr;
pr = nr ./ n;
pr(n == 0) = 0.5;

% build structure and return
cond_stats = struct('hs', hs, 'rt_corr', rt_corr, 'rt_incorr', rt_incorr, ...
    'rt_corr_var', rt_corr_var, 'rt_incorr_var', rt_incorr_var, ...
    'pr', pr, 'n', n, 'nr', nr, 'n_corr', n_corr, 'n_incorr', n_incorr);
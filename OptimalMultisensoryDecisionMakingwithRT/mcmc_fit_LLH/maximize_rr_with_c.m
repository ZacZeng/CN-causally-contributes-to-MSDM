function maximize_rr_with_c(subj_id, cost, restart_num, no_bias, model_name)
%% tunes bounds to maximize reward rate for given subject id, assuming cost
%
% The function computes the empirical reward rate and that predicted by the
% model and tunes the bound parameters of the latter to maximise this rate.
% Maximization is performed by assuming all model parameters other than the
% bounds to be fixed, and by performing gradient descent on the reward rate
% to find the bounds that maximize the latter. The function assumes a
% reward of 1 for correct choices, and a punishment of 0 for incorrect
% ones. cost is the cost for accumulating evidence, in units of rewards.
%
% The function allows for two different models (specified by model_name):
% - optim: assumes to find model fits for the optim model, which does not
%   assume any functional form on how the bounds depend on the coherence.
%   It uses the mf_optim fits as base parameters for the maximization. It
%   searches the best reward rate in the space of all possible bound
%   parameter combinations.
% - optim_powcomb: the same as for optim, only that the function uses the
%   optim_powcomb model, which has one bound parameter per condition over
%   which this function optimizes the bound parameters.
%
% For either model, the function runs gradient ascent with random restart a
% restart_num of times, using the parameter bounds of the corresponding
% model as initializer.
%
% If no_bias is given and true, that the biases of the model are set to 0 before
% computing the model reward rates.


%% detect if running on cluster, process arguments accordingly
on_cluster = ~usejava('desktop');
if nargin < 5, model_name = 'optim'; end
if on_cluster
    fprintf('Desktop not detected - assuming running on cluster\n');
    if nargin < 4, no_bias = false;
    else no_bias = logical(str2num(no_bias));
    end
    if nargin < 3, restart_num = 100;
    else restart_num = str2double(restart_num);
    end
    cost = str2double(cost);
    base_folder = '/home/jdrugowitsch/vis_vest/mcmc_fit_LLH';
else
    base_folder = '.';
    if nargin < 4, no_bias = false;
    else no_bias = logical(no_bias);
    end
    if nargin < 3, restart_num = 10; end
end
cost_str = sprintf('c%4.2f', cost);
if strcmp(model_name, 'optim'), model_fn = str2func('mf_optim');
elseif strcmp(model_name, 'optim_powcomb'), model_fn = str2func('mf_optim_powcomb');
else
    error('Unknown model name %s', model_name);
end
data_folder = [base_folder filesep 'fit_data'];
data_filename = [data_folder filesep 'mf_' model_name '_' subj_id '.mat'];
if no_bias
    rr_data_filename = [data_folder filesep 'max_rr_' cost_str '_no_bias_' model_name '_' subj_id '.mat'];
else
    rr_data_filename = [data_folder filesep 'max_rr_' cost_str '_' model_name '_' subj_id '.mat'];
end
opt=optimset('Display','iter','MaxFunEval',10000,'FunValCheck','on', ...
    'MaxIter',5000,'Algorithm','interior-point','TolX',1e-10);
delta_t = 5e-3;
t_max = 2;
R = [1 0];
iti = 6;
fprintf('Maximizing reward rate for subject %s, with %d random restarts, %s\n', ...
    subj_id, restart_num, cost_str);


%% load model fit data from optim model & compute various reward rates
dm = load(data_filename);
dstats = dm.dstats;
% reward rate for immediate random decisions
c = dm.cohs;
c_num = length(c);
if ~any(c_num == [3 6]), error('Unknown c_num detected'); end
p = dm.best_p;
if strcmp(model_name, 'optim')
    % optim model
    tnds = p((3*c_num+3):(3*c_num+5));
    % indices of bounds within parameter vector
    bp_idcs = [(c_num+1):(2*c_num) (2*c_num+2) (2*c_num+3):(3*c_num+2)];
    % bias_idx includes lapses
    if c_num == 3, bias_idx = 15:22;
    else bias_idx = 24:37;  end
else
    % optim_powcomb model
    tnds = p(10:12);
    bp_idcs = [3 6 8];
    if c_num == 3, bias_idx = 13:20;
    else bias_idx = 13:26;  end
end
cond_n = [0 sum(dstats.vest(1).n) 0];
for i = 1:length(dstats.vis), cond_n(1) = cond_n(1) + sum(dstats.vis(i).n); end
for i = 1:length(dstats.comb), cond_n(3) = cond_n(3) + sum(dstats.comb(i).n); end
mean_tnd = sum(tnds.*cond_n)/sum(cond_n);
rr_rnd = (0.5*R(1)+0.5*R(2)-mean_tnd*cost)/(mean_tnd+iti);
fprintf('random reward rate       : %6.4f\n\n', rr_rnd);
% empirical reward rate and untuned model reward rate
rr_emp = rrate(dstats, dstats, mean_tnd, R, cost, iti);
rr_ini = rrate(dm.mstats, dstats, mean_tnd, R, cost, iti);
fprintf('empirical reward rate    : %6.4f\n', rr_emp);
fprintf('initial model reward rate: %6.4f\n', rr_ini);


%% load subject data
if on_cluster
    ds = load([base_folder filesep 'subj_data' filesep subj_id '.mat']);
else
    ds = load_condi_data(subj_id, true);
end
n = dstats.n;
% add additional stats to d to be used by the model
ds.hs = dstats.vis(1).hs;
ds.delta_t = delta_t;
ds.t_max = t_max;
ds.subj_id = subj_id;
[ds.v_t, ds.a_t] = get_vel_acc_profile(delta_t);


%% extract previous bounds and find max/min values
pstruct = model_fn(ds);
p_min = pstruct.p_min;
p_max = pstruct.p_max;
p_names = pstruct.p_names;
p_num = length(p_min);
bp_min = p_min(bp_idcs);
bp_max = p_max(bp_idcs);
bp_ini = p(bp_idcs);
bp_names = p_names(bp_idcs);
bp_num = length(bp_idcs);


%% if no_bias option is chosen, re-compute initial model reward rate
if no_bias
    % bias_idx contains biases and lapse parameter
    p(bias_idx) = 0;
    rr_ini_no_bias = model_rrate(bp_ini, p, bp_idcs, ds, dstats, mean_tnd, R, cost, iti, model_fn);
    fprintf('model reward rate (no bias): %6.4f\n', rr_ini_no_bias);
end


%% reward-rate function
RR_func = @(bp) model_rrate(max(bp_min, min(bp_max, bp)),p,bp_idcs,ds,dstats,mean_tnd,R,cost,iti,model_fn);


%% find initial gradient / hessian
fprintf('\n----Estimating initial gradient and Hessian\n');
[grad_ini, grad_ini_err, grad_ini_final_delta] = gradest(RR_func, bp_ini);
fprintf('Found gradient  : %s\n', num2str(grad_ini));
fprintf('Estimation error: %s\n', num2str(grad_ini_err));
[H_ini, H_ini_err] = hessian(RR_func, bp_ini);
fprintf('Found Hessian   :\n');  disp(H_ini);
fprintf('Estimation error:\n');  disp(H_ini_err);


%% sample initial bound starting points based on previously best bound
% sampling is performed by first drawing from x~N(0,1) truncated to [-2,2],
% and then multiplying the associated bound parameter by exp(x log(10)/2).
% This causes the initial bound parameter to be randomly scaled by a scalar
% within the range [0.1, 10].
bp_inis = -3*ones(1, restart_num*bp_num);
bp_out_of_range = bp_inis < -2 | bp_inis > 2;
while any(bp_out_of_range)
    bp_inis(bp_out_of_range) = randn(1, sum(bp_out_of_range));
    bp_out_of_range = bp_inis < -2 | bp_inis > 2;
end
bp_inis = exp(reshape(bp_inis, restart_num, bp_num)*0.5*log(10));
bp_inis = bsxfun(@times,bp_inis,bp_ini);
bp_inis = bsxfun(@max, bp_min, bsxfun(@min, bp_max, bp_inis));


%% perform gradient ascent on reward rate from different starting points
rrs = zeros(1, restart_num);
ps = repmat(p,restart_num,1);
for i = 1:restart_num
    fprintf('\n----Gradient ascent #%d\n', i);
    [bp, rr] = fmincon(@(bp) -RR_func(bp),...
        bp_inis(i,:),[],[],[],[],bp_min,bp_max,[],opt);
    rr = -rr;
    ps(i,bp_idcs) = bp;
    rrs(i) = rr;
    fprintf('Found reward rate %6.4f\n', rr);
end
fprintf('Finished gradient descent runs\n\n');
% get model statistics for best reward rate
[rr_best,i] = max(rrs);
p_best = ps(i,:);
bp_best = p_best(bp_idcs);
fprintf('Best reward rate: %6.4f\n', rr_best);
[~,mstats] = RR_func(bp_best);
fprintf('at parameters\n');
for p_idx = 1:length(bp_idcs)
    fprintf('%s = %f\n', bp_names{p_idx}, p_best(bp_idcs(p_idx)));
end

%% find final gradient / hessian
fprintf('\n----Estimating final gradient and Hessian\n');
[grad_best, grad_best_err, grad_best_final_delta] = gradest(RR_func, bp_best);
fprintf('Found gradient  : %s\n', num2str(grad_best));
fprintf('Estimation error: %s\n', num2str(grad_best_err));
[H_best, H_best_err] = hessian(RR_func, bp_best);
fprintf('Found Hessian   :\n');  disp(H_best);
fprintf('Estimation error:\n');  disp(H_best_err);


%% writing data to file
fprintf('Writing data to %s\n', rr_data_filename);
if no_bias
    save(rr_data_filename, 'rr_emp', 'rr_ini', 'rr_rnd', 'rr_best', ...
        'rr_ini_no_bias', ...
        'c', 'R', 'iti', 'delta_t', 't_max', ...
        'bp_inis', 'rrs', 'ps', 'p_names', 'bp_idcs',...
        'dstats', 'mstats', 'subj_id', ...
        'grad_ini', 'grad_ini_err', 'grad_ini_final_delta', 'H_ini', 'H_ini_err', ...
        'grad_best', 'grad_best_err', 'grad_best_final_delta', 'H_best', 'H_best_err');
else
    save(rr_data_filename, 'rr_emp', 'rr_ini', 'rr_rnd', 'rr_best', ...
        'c', 'R', 'iti', 'delta_t', 't_max', ...
        'bp_inis', 'rrs', 'ps', 'p_names', 'bp_idcs',...
        'dstats', 'mstats', 'subj_id', ...
        'grad_ini', 'grad_ini_err', 'grad_ini_final_delta', 'H_ini', 'H_ini_err', ...
        'grad_best', 'grad_best_err', 'grad_best_final_delta', 'H_best', 'H_best_err');
end


%% exit matlab if running on cluster
if on_cluster
    exit;
end


function [rr, mstats] = model_rrate(bp, params, bp_idcs, ds, dstats, mean_tnd, R, cost, iti, model_fn)
%% returns the model reward rate for the given parameters
%
% bp are the bound parameters, params all parameters in combination, and
% bp_idcs determines the mapping of bp into params. ds is the subject data,
% dstats are some subject data statistics, R is the reward/punishment, and
% iti is the inter-trial interval.

% map bp into parameter vector
params(bp_idcs) = bp;
% compute model predictions and reward reward rate for these
mstats = model_fn(ds, params);
rr = rrate(mstats, dstats, mean_tnd, R, cost, iti);


function rr = rrate(mstats, dstats, mean_tnd, R, cost, iti)
%% returns the reward rate for given model stats and reward/ITI
%
% mstats is a structure containing the model statistics per condition.
% dstats is the same structure containing subjects data and more info on
% the task (like the number of repetitions per condition).
%
% R(1) is the reward given for correct decisions and R(2) the reward for
% incorrect decisions. iti is the inter-trial interval.

% compute mean RT and PC
mean_rt = 0;
mean_pc = 0;
n = 0;
% visual
for i = 1:length(mstats.vis)
    pc = mstats.vis(i).pr;
    pc(dstats.vis(i).hs < 0) = 1-pc(dstats.vis(i).hs < 0);
    rt = mstats.vis(i).rt_corr.*pc + mstats.vis(i).rt_incorr.*(1-pc);
    mean_rt = mean_rt + sum(rt.*dstats.vis(i).n);
    mean_pc = mean_pc + sum(pc.*dstats.vis(i).n);
    n = n + sum(dstats.vis(i).n);
end
% vestibular
pc = mstats.vest(1).pr;
pc(dstats.vest(1).hs < 0) = 1-pc(dstats.vest(1).hs < 0);
rt = mstats.vest(1).rt_corr.*pc + mstats.vest(1).rt_incorr.*(1-pc);
mean_rt = mean_rt + sum(rt.*dstats.vest(1).n);
mean_pc = mean_pc + sum(pc.*dstats.vest(1).n);
n = n + sum(dstats.vest(1).n);
% combined
for i = 1:length(mstats.comb)
    pc = mstats.comb(i).pr;
    pc(dstats.comb(i).hs < 0) = 1-pc(dstats.comb(i).hs < 0);
    rt = mstats.comb(i).rt_corr.*pc + mstats.comb(i).rt_incorr.*(1-pc);
    mean_rt = mean_rt + sum(rt.*dstats.comb(i).n);
    mean_pc = mean_pc + sum(pc.*dstats.comb(i).n);
    n = n + sum(dstats.comb(i).n);
end
% re-normalise
mean_rt = mean_rt / n;
mean_pc = mean_pc / n;

% reward rate is expected reward per trial over average trial time
rr = (mean_pc*R(1) + (1-mean_pc)*R(2) - max(0, mean_rt - mean_tnd) * cost) / (mean_rt+iti);


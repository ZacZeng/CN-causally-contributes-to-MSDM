function plot_rr_with_bound
%% plots how the reward rate changes with bound height

%% settings
% model and data files to use
model_fn = str2func('mf_optim');
fit_data_file = @(subj_id) ['fit_data' filesep 'mf_optim_' subj_id '.mat'];
rr_data_file = @(subj_id, cost) sprintf('fit_data%smax_rr_c%4.2f_optim_%s.mat', filesep, cost, subj_id);
% subjects to analyze + example subject for rr over bound plots.
subj_ids = {'231139','231140','548203','548214','950091','950107','967716','jason','eliana','kalpana'};
%subj_ids = {'231139','231140','548203'};
%ex_subj = '231139';
% evaluating rr over bound changes
bound_scaling = linspace(0.001, 3, 100);
% other settings
R = [1 0];
cs = [0 0.1 0.2];
ccols = [0 0 0; 0.4 0 0; 0.8 0 0];
iti = 6;
delta_t = 5e-3;
t_max = 2;

%% analyze parameters and rr per subject
subj_num = length(subj_ids);
best_bps = cell(subj_num, length(cs));
fitted_bps = cell(subj_num);
best_eucl_dist = zeros(subj_num, length(cs));
best_rr_dist = zeros(subj_num, length(cs));
emp_rrs = zeros(subj_num, length(cs));
best_rrs = zeros(subj_num, length(cs));
rnd_rrs = zeros(subj_num, length(cs));
rr_with_bound = zeros(subj_num, length(bound_scaling), length(cs));
fprintf('Processing %d subjs\n', subj_num);
for subj_idx = 1:subj_num
    subj_id = subj_ids{subj_idx};
    fprintf('%s ... ', subj_id);
    % load subject fits
    dm = load(fit_data_file(subj_id));
    dstats = dm.dstats;
    best_p = dm.best_p;
    % load subject data stats
    ds = load_condi_data(subj_id, false);
    % add additional stats to d to be used by the model
    ds.hs = dstats.vis(1).hs;
    ds.delta_t = delta_t;
    ds.t_max = t_max;
    ds.subj_id = subj_id;
    [ds.v_t, ds.a_t] = get_vel_acc_profile(delta_t);
    % average non-decision time
    coh_num = length(dm.cohs);
    if ~any(coh_num == [3 6]), error('Unknown coh_num detected'); end
    % optim model - THE NEXT TWO LINES BREAK IF ANOTHER MODEL THAN optim IS USED
    tnds = best_p((3*coh_num+3):(3*coh_num+5));
    % indices of bounds within parameter vector
    bp_idcs = [(coh_num+1):(2*coh_num) (2*coh_num+2) (2*coh_num+3):(3*coh_num+2)];
    % bias_idx includes lapses
    cond_n = [0 sum(dstats.vest(1).n) 0];
    for i = 1:length(dstats.vis), cond_n(1) = cond_n(1) + sum(dstats.vis(i).n); end
    for i = 1:length(dstats.comb), cond_n(3) = cond_n(3) + sum(dstats.comb(i).n); end
    mean_tnd = sum(tnds.*cond_n)/sum(cond_n);
    % extract fitted bound parameters
    fitted_bps{subj_idx} = best_p(bp_idcs);
    % extract min/max parameter values, imposed by model
    pstruct = model_fn(ds);
    bp_min = pstruct.p_min(bp_idcs);
    bp_max = pstruct.p_max(bp_idcs);
    
    for c_idx = 1:length(cs)
        c = cs(c_idx);
        fprintf('%4.2f ', c);
        % extract highest-rr parameters
        d = load(rr_data_file(subj_id, c));
        assert(all(bp_idcs == d.bp_idcs));
        best_rrs(subj_idx, c_idx) = d.rr_best;
        rnd_rrs(subj_idx, c_idx) = d.rr_rnd;
        emp_rrs(subj_idx, c_idx) = d.rr_emp;
        [~,i] = max(d.rrs);
        best_bps{subj_idx, c_idx} = d.ps(i,bp_idcs);

        % eucledian distance to optimal bound (<1 - below, >1 - above
        best_rel = bound_to_rel(best_bps{subj_idx, c_idx});
        fitted_rel = bound_to_rel(fitted_bps{subj_idx});
        best_eucl_dist(subj_idx, c_idx) = (best_rel * fitted_rel') / (best_rel * best_rel');
        % compute rr over bound, in rel space
        %emp_rrs(subj_idx, c_idx) = rrate(dstats, dstats, R, c, iti);
        for i = 1:length(bound_scaling)
            bps = rel_to_bound(bound_scaling(i) * best_rel);
            rr_with_bound(subj_idx, i, c_idx) = ...
                model_rrate(max(bp_min, min(bp_max, bps)), ...
                            best_p, bp_idcs, ds, dstats, mean_tnd, R, c, iti, model_fn);
        end
        % find point on rr curve at which overlap
        if emp_rrs(subj_idx, c_idx) > max(rr_with_bound(subj_idx, :, c_idx))
            % empirical rr exceeds best one found for model -> set empirial
            % point on maximum rr found by model
            [emp_rrs(subj_idx, c_idx),i] = max(rr_with_bound(subj_idx, :, c_idx));
            best_rr_dist(subj_idx, c_idx) = bound_scaling(i);
        else
            if best_eucl_dist(subj_idx) <= 1
                % below 'optimal' settings
                max_idx = find(bound_scaling <= 1, 1, 'last');
                best_rr_dist(subj_idx, c_idx) = interp1(rr_with_bound(subj_idx, 1:max_idx, c_idx), ...
                    bound_scaling(1:max_idx), emp_rrs(subj_idx, c_idx), 'cubic');
            else
                % above 'optimal' settings
                min_idx = find(bound_scaling > 1, 1, 'first');
                best_rr_dist(subj_idx, c_idx) = interp1(rr_with_bound(subj_idx, min_idx:end, c_idx), ...
                    bound_scaling(min_idx:end), emp_rrs(subj_idx, c_idx), 'cubic', 'extrapolate');
            end
        end
    end
    fprintf('\n');
end

%% plot distance to bound per subject
pbar = 4/3;
figure('Color', 'white'); hold on;
max_dist = max(max(best_rr_dist));
if max_dist < 1.2; xlim([0 1.2]); end
ylim([0.5 (subj_num+0.5)]);
plot([1 1], ylim, 'k--');
for c_idx = 1:length(cs);
    plot(best_rr_dist(:,c_idx), 1:subj_num, 'o', ...
        'MarkerFaceColor', ccols(c_idx, :), 'MarkerEdgeColor', 'none', 'MarkerSize', 4);
end
xlabel('subject bound compared to optimal');
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
    'YTick',[]);
%plot(1:subj_num, best_eucl_dist, 'ko');
%plot(1:subj_num, best_rr_dist(:,1), 'ro');

pbar = 4/3;
figure('Color', 'white'); hold on;
max_dist = max(max(best_eucl_dist));
if max_dist < 1.2; xlim([0 1.2]); end
ylim([0.5 (subj_num+0.5)]);
plot([1 1], ylim, 'k--');
for c_idx = 1:length(cs);
    plot(best_eucl_dist(:,c_idx), 1:subj_num, 'o', ...
        'MarkerFaceColor', ccols(c_idx, :), 'MarkerEdgeColor', 'none', 'MarkerSize', 4);
end
xlabel('subject bound compared to optimal (eucledian)');
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
    'YTick',[]);


% plot real parameter per subject
for subj_idx = 1:subj_num
    pbar = 4/3;
    figure('Color','white'); hold on;  xlim([0 max(bound_scaling)]);
    for c_idx = 1:length(cs)
        plot(bound_scaling, rr_with_bound(subj_idx, :, c_idx), 'Color', ccols(c_idx, :));
        plot(best_rr_dist(subj_idx, c_idx), emp_rrs(subj_idx, c_idx), 'o', ...
            'MarkerFaceColor', ccols(c_idx, :), 'MarkerEdgeColor', 'none', 'MarkerSize', 4);
        plot(1, max(rr_with_bound(subj_idx, :, c_idx)), 'x', ...
            'Color', ccols(c_idx, :), 'MarkerSize', 4);
        set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
            'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
    end
end




function x = bound_to_rel(bound)
%% mapping bound into relevant variable
x = bound;

function bound = rel_to_bound(x)
%% mapping relevant variable to bound
bound = x;


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


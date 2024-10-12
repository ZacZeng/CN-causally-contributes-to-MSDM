function plot_rr_gradient
%% performs analysis of the reward rate gradient at the estimated bound

%% settings
fit_data_file = @(subj_id) ['fit_data' filesep 'mf_optim_' subj_id '.mat'];
rr_data_file = @(subj_id, cost) sprintf('fit_data%smax_rr_c%4.2f_optim_%s.mat', filesep, cost, subj_id);
subj_ids = {'231139','231140','548203','548214','950091','950107','967716','jason','eliana','kalpana'};
costs = [0.0 0.1 0.2];
subj_num = length(subj_ids);
cost_num = length(costs);
msize = 4;


%% load data
bound_diff = cell(subj_num, cost_num);
bound_dist = cell(subj_num, cost_num);
RR_curvature = cell(subj_num, cost_num);
RR_H = cell(subj_num, cost_num);
RR_grad = cell(subj_num, cost_num);
RR_H_eig = cell(subj_num, cost_num);
RR_H_dist = cell(subj_num, cost_num);
RR_H_cols = cell(subj_num, cost_num);
RR_obs_res = cell(subj_num, cost_num);
RR_pred_res = cell(subj_num, cost_num);
for subj_idx = 1:subj_num
    subj_id = subj_ids{subj_idx};
    for cost_idx = 1:cost_num
        % hessian / gradient at estimated bound location
        d = load(rr_data_file(subj_id, costs(cost_idx)));
        bp_idcs = d.bp_idcs;
        [~,i] = max(d.rrs);
        bps_best = d.ps(i,bp_idcs);
        RR_H{subj_idx, cost_idx} = d.H_ini;
        RR_grad{subj_idx, cost_idx} = d.grad_ini;
        RR_curvature{subj_idx, cost_idx} = diag(d.H_ini);
        RR_obs_res{subj_idx, cost_idx} = max(d.rrs) - d.rr_emp;
        % parameters for estimated bound location
        d = load(fit_data_file(subj_id));
        bps_fit = d.best_p(bp_idcs);
        bound_diff{subj_idx, cost_idx} = bps_best - bps_fit;
        bound_dist{subj_idx, cost_idx} = abs(bps_best - bps_fit);
        % project distance into Hessian subspace
        [V,D] = eig(RR_H{subj_idx, cost_idx});
        RR_H_eig{subj_idx, cost_idx} = diag(D)';
        % project distance onto eigenvectors (does not change norm)
        proj_dist = V' * (bps_best - bps_fit)';
        RR_H_dist{subj_idx, cost_idx} = proj_dist;
        % projected colors (and fix signs)
        cols = V' * getcolmat(size(V,1));
        for n = 1:size(cols,1)
            [~,i] = max(abs(cols(n,:)));
            if cols(n,i) < 0; cols(n,:) = -cols(n,:); end
            cols(n,:) = min(1.0, max(0.0, cols(n,:)));
        end
        RR_H_cols{subj_idx, cost_idx} = cols;
        % compute predicted reward rate residual, using 2nd order Taylor
        % expansion
        RR_pred_res{subj_idx, cost_idx} = RR_obs_res{subj_idx, cost_idx} + 0.0005 * (...
            - bound_diff{subj_idx, cost_idx} * RR_grad{subj_idx, cost_idx}' ...
           - 0.5 * ...
           bound_diff{subj_idx, cost_idx} * RR_H{subj_idx, cost_idx} * ...
           bound_diff{subj_idx, cost_idx}');
    end
end


%% plot results
pbar = 4/3;
% bound distance vs. curvature
for cost_idx = 1:cost_num
    fprintf('Cost %4.2f\n', costs(cost_idx));
    figure('Color', 'white');  hold on;
    all_dist = [];
    all_curvature = [];
    for subj_idx = 1:subj_num
        cond_num = length(bound_dist{subj_idx, cost_idx});
        colmat = getcolmat(cond_num);
        for i = 1:cond_num
            plot(bound_dist{subj_idx, cost_idx}(i), RR_H{subj_idx, cost_idx}(i,i), 'o', ...
                'MarkerFaceColor', colmat(i,:), 'MarkerEdgeColor', 'none', ...
                'MarkerSize', msize);
        end        
        %plot(bound_dist{subj_idx, cost_idx}, diag(RR_H{subj_idx, cost_idx}), 'o', ...
        %    'MarkerFaceColor', colorget(subj_idx), 'MarkerEdgeColor', 'none');
        all_dist = cat(1, all_dist, bound_dist{subj_idx, cost_idx}');
        all_curvature = cat(1, all_curvature, diag(RR_H{subj_idx, cost_idx}));
    end
    xlabel('absolute bound distance');
    ylabel(sprintf('curvature at estimated bound location, c=%4.2f',costs(cost_idx)));
    [rho,pval] = corr(all_dist, all_curvature, 'type', 'Spearman');
    fprintf('  rho = %6.4f, p = %6.4f\n', rho, pval);
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
       'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
end
% gradient vs. curvature
for cost_idx = 1:cost_num
    figure('Color', 'white');  hold on;
    for subj_idx = 1:subj_num
        cond_num = length(RR_grad{subj_idx, cost_idx});
        colmat = getcolmat(cond_num);
        for i = 1:cond_num
            plot(RR_grad{subj_idx, cost_idx}(i), RR_H{subj_idx, cost_idx}(i,i), ...
                'o', 'MarkerFaceColor', colmat(i,:), 'MarkerEdgeColor', 'none', ...
                'MarkerSize', msize);
        end
    end
    xlabel('gradient');
    ylabel(sprintf('curvature at estimated bound location, c=%4.2f',costs(cost_idx)));
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
      'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
end
% gradient vs. distance
for cost_idx = 1:cost_num
    figure('Color', 'white');  hold on;
    for subj_idx = 1:subj_num
        cond_num = length(bound_dist{subj_idx, cost_idx});
        colmat = getcolmat(cond_num);
        for i = 1:cond_num
            plot(bound_dist{subj_idx, cost_idx}(i), abs(RR_grad{subj_idx, cost_idx}(i)), ...
                'o', 'MarkerFaceColor', colmat(i,:), 'MarkerEdgeColor', 'none', ...
                'MarkerSize', msize);
        end
    end
    xlabel('absolute bound distance');
    ylabel(sprintf('gradient at estimated bound location, c=%4.2f',costs(cost_idx)));
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
      'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
end
% bound distance vs. curvature, rotated space
for cost_idx = 1:cost_num
    figure('Color', 'white');  hold on;
    for subj_idx = 1:subj_num
        cond_num = length(RR_H_dist{subj_idx, cost_idx});
        colmat = RR_H_cols{subj_idx, cost_idx};
        for i = 1:cond_num
            plot(abs(RR_H_dist{subj_idx, cost_idx}(i)), RR_H_eig{subj_idx, cost_idx}(i), ...
                'o', 'MarkerFaceColor', colmat(i,:), 'MarkerEdgeColor', 'none', ...
                'MarkerSize', msize);
        end
    end
    xlabel('projected bound distance');
    ylabel(sprintf('curvature at estimated bound location, c=%4.2f',costs(cost_idx)));
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
      'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
end

%% residual bounds, observed vs. predicted, cost 0.0
pbar = 1;
figure('Color', 'white');  hold on;
xlim([0 0.014]);  ylim([0 0.014]);
plot(xlim, ylim, '-', 'LineWidth', 0.5, 'Color', [1 1 1]*0.5);
plot([RR_obs_res{:,1}], [RR_pred_res{:,1}], 'ko', ...
    'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 5);
xlabel('observed optimal - measured reward rate');
ylabel('predicted optimal - measured reward rate');
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
      'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));


%% plot colors
colmat = getcolmat(13);
figure('Color', 'white');  hold on;
for n = 1:6
    plot(1, n, 'o', ...
        'MarkerFaceColor', colmat(n,:), 'MarkerEdgeColor', 'none', 'MarkerSize', msize);
    plot(3, n, 'o', ...
        'MarkerFaceColor', colmat(n+7,:), 'MarkerEdgeColor', 'none', 'MarkerSize', msize);
end
plot(2, 1, 'o', ...
     'MarkerFaceColor', colmat(7,:), 'MarkerEdgeColor', 'none', 'MarkerSize', msize);


function colmat = getcolmat(cond_num)
%% returns the color matrix for the given number of conditions
%
% the matrix is in [vis; vest; comb] order. The cond_num gives the total
% number of conditions

if cond_num == 7
    colmat = [0  50 200; 0 100 200; 0 180 200; ...
              0 200   0; ...
              200  50 0; 200 100 0; 200 180 0] / 255;
elseif cond_num == 13
    colmat = [0   0 200; 0  25 200; 0  50 200; ...
              0 100 200; 0 140 200; 0 180 200;
              0 200   0; ...
              200   0 0; 200  25 0; 200  50 0; ...
              200 100 0; 200 140 0; 200 180 0] / 255;
else
    error('Unknown number of conditions');
end

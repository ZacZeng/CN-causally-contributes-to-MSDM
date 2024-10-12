function plot_sample_stats(subj_id, model_id, burnin)
%% plots some statistics of the samples for the given model

if nargin < 3, burnin = []; end
conf_lo = 2.5;
conf_up = 97.5;

%% load samples & remove burn-in
d = load(['fit_data' filesep model_id '_' subj_id '.mat']);
if isempty(burnin), burnin = floor(length(d.llhs)/2); end
d.s = d.s((burnin+1):end,:);
d.llhs = d.llhs((burnin+1):end,:);
trials = length(d.llhs);
hist_bins = min(200, floor(trials/20));
p_num = length(d.p_names);
[~,s_mode] = max(d.llhs);


%% distribution of llhs
figure;
hist(d.llhs, hist_bins);
hold on;
plot(ones(1,2)*d.best_llh, get(gca,'YLim'),'k-');
plot(ones(1,2)*max(d.llhs), get(gca,'YLim'),'k--');
xlabel('log-likelihood');


%% distribution of parameters
for i = 1:p_num
    figure;
    p_lo = prctile(d.s(:,i),conf_lo);
    p_up = prctile(d.s(:,i),conf_up);
    p_mode = d.best_p(i);
    p_mode_s = d.s(s_mode,i);
    p_mean_s = mean(d.s(:,i));
    p_sd = sqrt(var(d.s(:,i)));
    p_lo_s = p_mean_s - 1.96 * p_sd;
    p_up_s = p_mean_s + 1.96 * p_sd;
    
    % determine plotting range
    min_lo = min([p_lo p_mode p_mode_s p_mean_s p_lo_s]);
    max_up = max([p_up p_mode p_mode_s p_mean_s p_up_s]);
    p_span = max_up - min_lo;
    min_lo = min_lo - p_span/4;
    max_up = max_up + p_span/4;
    
    % plot histogram
    subplot(2,1,1);  hold on;
    hist(d.s(:,i), hist_bins);
    ylims = get(gca,'YLim');
    if p_mode < max(p_lo, p_lo_s) || p_mode > min(p_up, p_up_s)
        plot(ones(1,2)*p_mode,ylims,'k-','LineWidth',2);
    else
        plot(ones(1,2)*p_mode,ylims,'k-');
    end
    plot(ones(1,2)*p_lo,ylims,'r-');
    plot(ones(1,2)*p_up,ylims,'g-');
    if p_mode_s < max(p_lo, p_lo_s) || p_mode_s > min(p_up, p_up_s)
        plot(ones(1,2)*p_mode_s,ylims,'k--','LineWidth',2);
    else
        plot(ones(1,2)*p_mode_s,ylims,'k--');
    end
    if p_mean_s < max(p_lo, p_lo_s) || p_mean_s > min(p_up, p_up_s)
        plot(ones(1,2)*p_mean_s,ylims,'k:','LineWidth',2);
    else
        plot(ones(1,2)*p_mean_s,ylims,'k:');
    end
    plot(ones(1,2)*p_lo_s,ylims,'r--');
    plot(ones(1,2)*p_up_s,ylims,'g--');
    xlim([min_lo max_up]);
    title(d.p_names{i});
    
    % plot samples
    subplot(2,1,2);  hold on;
    plot(d.s(:,i),1:trials,'b-');
    ylims = [1 trials];
    if p_mode < max(p_lo, p_lo_s) || p_mode > min(p_up, p_up_s)
        plot(ones(1,2)*p_mode,ylims,'k-','LineWidth',2);
    else
        plot(ones(1,2)*p_mode,ylims,'k-');
    end
    plot(ones(1,2)*p_lo,ylims,'r-');
    plot(ones(1,2)*p_up,ylims,'g-');
    if p_mode_s < max(p_lo, p_lo_s) || p_mode_s > min(p_up, p_up_s)
        plot(ones(1,2)*p_mode_s,ylims,'k--','LineWidth',2);
    else
        plot(ones(1,2)*p_mode_s,ylims,'k--');
    end
    if p_mean_s < max(p_lo, p_lo_s) || p_mean_s > min(p_up, p_up_s)
        plot(ones(1,2)*p_mean_s,ylims,'k:','LineWidth',2);
    else
        plot(ones(1,2)*p_mean_s,ylims,'k:');
    end
    plot(ones(1,2)*p_lo_s,ylims,'r--');
    plot(ones(1,2)*p_up_s,ylims,'g--');
    xlim([min_lo max_up]);
    ylim(ylims);
end

function plot_rr
%% plot reward rate fits/info

%% settings
% one bound per coherence
%data_file = @(subj_id) ['fit_data' filesep 'max_rr_optim_' subj_id '.mat'];
data_file = @(subj_id, cost) sprintf('fit_data%smax_rr_c%4.2f_optim_%s.mat', filesep, cost, subj_id);
costs = [0.0 0.1 0.2];
data_file_nb = @(subj_id) ['fit_data' filesep 'max_rr_no_bias_optim_', subj_id '.mat'];
% using bound function (restricted)
data_file_res = @(subj_id) ['fit_data' filesep 'max_rr_optim_powcomb_' subj_id '.mat'];
data_file_res_nb = @(subj_id) ['fit_data' filesep 'max_rr_no_bias_optim_powcomb_', subj_id '.mat'];
subj_ids = {'231139','231140','548203','548214','950091','950107','967716','jason','eliana','kalpana'};
model_names = {'c0.0', 'c0.1', 'c0.2', 'constr', 'nb', 'constr nb'};
base_model_idx = 1;
comp_model_idx = 4;

%% collect data
subj_num = length(subj_ids);
cost_num = length(costs);
rr_emp = NaN(cost_num,subj_num);
rr_rnd = NaN(3 + cost_num,subj_num);
rr_best = NaN(3 + cost_num,subj_num);
for subj_idx = 1:subj_num
    % 1..cost_num = one bound per coherence, biased
    for c_idx = 1:cost_num
        d = load(data_file(subj_ids{subj_idx},costs(c_idx)));
        rr_emp(c_idx,subj_idx) = d.rr_emp;
        rr_rnd(c_idx,subj_idx) = d.rr_rnd;
        rr_best(c_idx,subj_idx) = d.rr_best;
    end
    % cost_num + 1 = restricted bounds, biased
    d = load(data_file_res(subj_ids{subj_idx}));
    rr_rnd(cost_num+1,subj_idx) = d.rr_rnd;
    rr_best(cost_num+1,subj_idx) = d.rr_best;
    % cost_num + 2 = one bound per coherence, non-biased
    d = load(data_file_nb(subj_ids{subj_idx}));
    rr_rnd(cost_num+2,subj_idx) = d.rr_rnd;
    rr_best(cost_num+2,subj_idx) = d.rr_best;
    % cost_num + 3 = restricted bounds, non-biased
    d = load(data_file_res_nb(subj_ids{subj_idx}));
    rr_rnd(cost_num+3,subj_idx) = d.rr_rnd;
    rr_best(cost_num+3,subj_idx) = d.rr_best;
end
rel_rr_emp = NaN(3 + cost_num, subj_num);
rel_rr_emp(1:cost_num, :) = rr_emp ./ rr_best(1:cost_num, :);
rel_rr_emp((cost_num+1):end, :) = bsxfun(@rdivide, rr_emp(1,:), rr_best((cost_num+1):end, :));
%rel_rr_emp = bsxfun(@rdivide, rr_emp, rr_best);
rel_rr_rnd = rr_rnd ./ rr_best;


%% plot data
pbar = 4/3;
figure('Color','white');  hold on;
xlim([0.5 (subj_num+0.5)]);  ylim([0 1]);
plot(xlim,[1 1]*0.9,'-','Color',[1 1 1]*0.5);
plot(xlim,[1 1]*0.95,'-','Color',[1 1 1]*0.5);
for subj_idx = 1:subj_num
    rectangle('Position',[(subj_idx-0.4) 0 0.375 rel_rr_emp(base_model_idx,subj_idx)],...
        'EdgeColor','none','FaceColor',[0 0 200]/255);
    rectangle('Position',[(subj_idx+0.0025) 0 0.375 rel_rr_rnd(base_model_idx,subj_idx)],...
        'EdgeColor','none','FaceColor',[200 0 0]/255);
end
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
    'XTick',[],...
    'YTick',0:0.1:1,'YTickLabel',{'0','','0.2','','0.4','','0.6','','0.8','','1'});
xlabel('subject');
ylabel('relative reward rate');


%% plot data, no bias
pbar = 4/3;
figure('Color','white');  hold on;
xlim([0.5 (subj_num+0.5)]);  ylim([0 1]);
plot(xlim,[1 1]*0.9,'-','Color',[1 1 1]*0.5);
plot(xlim,[1 1]*0.95,'-','Color',[1 1 1]*0.5);
for subj_idx = 1:subj_num
    rectangle('Position',[(subj_idx-0.4) 0 0.375 rel_rr_emp(comp_model_idx,subj_idx)],...
        'EdgeColor','none','FaceColor',[0 0 200]/255);
    rectangle('Position',[(subj_idx+0.0025) 0 0.375 rel_rr_rnd(comp_model_idx,subj_idx)],...
        'EdgeColor','none','FaceColor',[200 0 0]/255);
end
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
    'XTick',[],...
    'YTick',0:0.1:1,'YTickLabel',{'0','','0.2','','0.4','','0.6','','0.8','','1'});
xlabel('subject');
ylabel('relative reward rate (no bias)');


%% compare optimal rate, no bias vs. bias
pbar = 1;
figure('Color','white');  hold on;
xlim([0.11 0.13]);  ylim([0.11 0.13]);
plot(xlim, ylim, 'k--');
plot(rr_best(base_model_idx,:), rr_best(comp_model_idx,:), 'ko', 'LineWidth', 1, 'MarkerFaceColor', [1 1 1],...
    'MarkerSize', 6);
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
    'XTick', [0.11 0.12 0.13], 'YTick', [0.11 0.12 0.13]);
xlabel('best reward rate, bias');
ylabel('best reward rate, no bias');

pbar = 4/3;
figure('Color','white');  hold on;
xlim([0.5 (subj_num+0.5)]);  ylim([0 1]);
for subj_idx = 1:subj_num
    rectangle('Position',[(subj_idx-0.4) 0 0.8 (rr_best(base_model_idx,subj_idx)/rr_best(comp_model_idx,subj_idx))],...
        'EdgeColor','none','FaceColor',0.5*[1 1 1]);
end
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
    'XTick',[],...
    'YTick',0:0.1:1,'YTickLabel',{'0','','0.2','','0.4','','0.6','','0.8','','1'});
ylabel('best RR (bias) / best RR (no bias)');
xlabel('subjects');


%% compare different optimality criteria with box plot
pbar = 4/3;
figure('Color','white');  hold on; xlim([0.5 (length(model_names)+0.5)]);  ylim([0 1]);
rel_rr_emp_median = prctile(rel_rr_emp(1,:),50);
rel_rr_rnd_median = prctile(rel_rr_rnd(1,:),50);
plot(xlim, [1 1]*rel_rr_emp_median, '--', 'Color', [0 0 0.8]);
plot(xlim, [1 1]*rel_rr_rnd_median, '--', 'Color', [0.8 0 0]);
for m_idx = 1:length(model_names)
    plot_box(m_idx-0.15,0.28,0.1,rel_rr_emp(m_idx,:),[0.4 0.4 0.9],[0 0 0.8]);
    plot_box(m_idx+0.15,0.28,0.1,rel_rr_rnd(m_idx,:),[0.9 0.4 0.4],[0.8 0 0]);
    %[~,p,~,stats] = ttest(rel_rr_emp(m_idx,:), rel_rr_rnd(m_idx,:));
    %fprintf('%10s, t(%d) = %f, p = %f\n', model_names{m_idx}, stats.df, stats.tstat, p);
    [p,~,stats] = signrank(rel_rr_emp(m_idx,:), rel_rr_rnd(m_idx,:));
    fprintf('%10s, srank = %f, p = %f\n', model_names{m_idx}, stats.signedrank, p);
end
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
    'XTick',1:length(model_names),'XTickLabel',model_names,...
    'YTick',0:0.1:1,'YTickLabel',{'0','','0.2','','0.4','','0.6','','0.8','','1'});
ylabel('relative reward rate');


function plot_box(x,wb,wt,d,cb,cl)
%% creates a box-plot at x of data d (vector).
% wb is the with of the rectangle (25-75 pctile), and wt the width of the
% top and bottom bars (showing minimum and maximum). cb is the overall
% color, and cl is color of median, both as [r g b].

ymax = max(d);
ymin = min(d);
y25 = prctile(d,25);
y50 = prctile(d,50);
y75 = prctile(d,75);
rectangle('Position',[(x-wb/2) y25 wb (y75-y25)],'EdgeColor','none','FaceColor',cb);
plot([1 1]*x,[ymin y25],'Color',cb);
plot([1 1]*x,[y75 ymax],'Color',cb);
plot(x+[-0.5 0.5]*wt,[1 1]*ymax,'Color',cb);
plot(x+[-0.5 0.5]*wt,[1 1]*ymin,'Color',cb);
plot(x+[-0.5 0.5]*wb,[1 1]*y50,'Color',cl);
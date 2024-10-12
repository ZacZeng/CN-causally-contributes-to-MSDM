function plot_fit_measure
%% plots the fit quality per subject


%% settings
data_folder = 'fit_data';
data_filename = @(subj_id, model_id)...
    [data_folder filesep model_id '_' subj_id '.mat'];


%% set of models and subjects
% model comparison, powcomb bounds
model_ids = {'mf_optim_powcomb','mf_sepk_powcomb',...
    'mf_nocuew_powcomb','mf_accw_powcomb','mf_velow_powcomb',...
    'mf_notempw_powcomb','mf_nocuetempw_powcomb'};
model_short = {'opt','sep','nc','acc','velo','nw','nw_nc'};

% same as above, without velocity-weighted model
%model_ids = {'mf_optim_powcomb','mf_sepk_powcomb',...
%    'mf_nocuew_powcomb','mf_accw_powcomb',...
%    'mf_notempw_powcomb','mf_nocuetempw_powcomb'};
%model_short = {'opt','sep','nc','acc','nw','nw_nc'};

% powcomb bound vs. variants (no bound, powsep, ...)
%model_ids = {'mf_optim_powcomb','mf_optim_powsep','mf_optim',...
%   'mf_optim_powcomb_onetnd','mf_optim_powcomb_nobias','mf_optim_powcomb_nobias_onetnd'};
%model_short = {'opt','sep','nb','onetnd','nobias','onetnd.nobias'};

% model comparison, powsep bounds
%model_ids = {'mf_optim_powsep','mf_sepk_powsep',...
%    'mf_nocuew_powsep','mf_notempw_powsep','mf_nocuetempw_powsep',...
%    'mf_accw_powsep','mf_velow_powsep'};
%model_short = {'opt','sep','nc','nw','nw_nc','acc','velo'};

% model comparison , no parametric bounds
%model_ids = {'mf_optim','mf_sepk','mf_nocuew','mf_notempw','mf_nocuetempw','mf_accw','mf_velow'};
%model_short = {'opt','sep','nc','nw','nw_nc','acc','velo'};


subj_ids = {'231139','231140','548203','548214','950091','950107','967716','sub01','subj02','subj03'};

subj_colors = [255 0 0; 0  255 0; 0 0 255; 255 0 255; 0 255  255; ...
    255 255 0; 0 0 0; 112 219 147; 181 166 66; 95 159 159; ...
    184 115 51; 47 79 47; 153 50 205; 135 31 120; 133 94 66; ...
    84 84 84; 142 35 35; 245 204 176; 35 142 35; 205 127 50; ...
    219 219 112; 192 192 192; 82 127 118; 159 159 95; 142 35 107; ...
    47 47 79; 235 199 158; 207 181 59; 255 127 0; 219 112 219; ...
    217 217 243; 89 89 171; 140 23 23; 35 142 104; 107 66 38]*(200/255/255);


%% collect fit quality
model_num = length(model_ids);
subj_num = length(subj_ids);
fit_llh = zeros(subj_num, model_num);
fit_bic = zeros(subj_num, model_num);
fit_r2 = zeros(subj_num, model_num);
fit_r2adj = zeros(subj_num, model_num);
for model_idx = 1:model_num
    model_id = model_ids{model_idx};
    model_fn = str2func(model_id);
    for subj_idx = 1:subj_num
        subj_id = subj_ids{subj_idx};
        d = load(data_filename(subj_id, model_id));
        fit_llh(subj_idx, model_idx) = d.best_llh;
        fit_bic(subj_idx, model_idx) = d.best_llh-...
            0.5*length(d.p_names)*log(d.dstats.n);
        fit_r2(subj_idx, model_idx) = get_R2(subj_id, model_fn);
        fit_r2adj(subj_idx, model_idx) = 1-(1-fit_r2(subj_idx,model_idx))*...
            (d.dstats.n-1)/(d.dstats.n-length(d.p_names)-1);
    end
    % output model quality
    fprintf('LLH_%s = c(%s)\n', model_short{model_idx}, ...
        cell2mat(...
            cat(2,arrayfun(@(i) sprintf('%8.6f, ',fit_llh(i,model_idx)),...
                           1:(subj_num-1),'UniformOutput',false),...
                  {sprintf('%8.6f',fit_llh(subj_num,model_idx))})));
end
% output code to create data frame for R
fprintf('LLH = c(%s)\n', ...
    cell2mat(cat(2,arrayfun(@(i) ['LLH_' model_short{i} ', '],...
                            1:(model_num-1),'UniformOutput',false),...
                   ['LLH_' model_short{length(model_short)}])));
fprintf('subj <- factor(rep(c(%s), %d))\n', ...
    cell2mat(cat(2,arrayfun(@(i) ['"' subj_ids{i} '", '],...
                            1:(subj_num-1),'UniformOutput',false),...
                   ['"' subj_ids{subj_num} '"'])), model_num);
fprintf('cond <- factor(rep(c(%s), each=%d))\n', ...
    cell2mat(cat(2,arrayfun(@(i) ['"' char(i-1+'a') '_' model_short{i} '", '],...
                            1:(model_num-1),'UniformOutput',false),...
                   ['"' char(model_num-1+'a') '_' model_short{model_num} '"'])),...
    subj_num);
fprintf('LLHd <- data.frame(LLH = LLH, subj=subj, cond=cond)\n');
% compute relative fits
rel_fit_llh = bsxfun(@minus, fit_llh, fit_llh(:,1));
rel_fit_bic = bsxfun(@minus, fit_bic, fit_bic(:,1));

% perform mixed-effects comparison
fprintf('\n\nPerforming mixed-effects comparison\n');
fam_num = model_num;
fam_model_idx = arrayfun(@(i) i, 1:fam_num, 'UniformOutput', false);
f_samp = family_wise(fit_bic, fam_model_idx, true, 20000);
mean_fam_p = mean(f_samp, 1);
sd_fam_p = sqrt(var(f_samp, 1));
[~, best_fam_idx] = max(mean_fam_p);
best_samp = true(size(f_samp,1),1);
for fam_idx = 1:fam_num
    if fam_idx == best_fam_idx, continue; end
    best_samp = best_samp & (f_samp(:,best_fam_idx) > f_samp(:,fam_idx));
end
exc_prob = mean(best_samp);
fprintf('Winning family: %s with exceedance probability %8.6f\n', ...
    model_ids{best_fam_idx}, exc_prob);
% plot family-wise comparison
pbar = 1;
figure('Color','white');  hold on;  xlim([0 1]);  ylim([0.5 (fam_num+0.5)]);
barh(mean_fam_p, 'FaceColor',[1 1 1]*0.5,'EdgeColor','none');
for fam_idx = 1:fam_num
    plot(mean_fam_p(fam_idx)+[-1 1]*sd_fam_p(fam_idx), [1 1]*fam_idx, 'k-');
end
ylabel('p(model)');
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
    'YTick',1:fam_num,'YTickLabel',model_short,'YDir','reverse');



%% basic plots
% individual fit likelihoods
figure;
bar(fit_llh);
set(gca,'XTick',1:subj_num,'XTickLabel',subj_ids);
xlabel('R2');

% individual relative BIC measures
figure;
bar(rel_fit_bic);
set(gca,'XTick',1:subj_num,'XTickLabel',subj_ids);
xlabel('LLH - LLH optim');
legend(model_short);

% individual relative BIC measures, again per subject
figure;
subplot(2,1,1);
bar(rel_fit_bic(:,2:end)'/log(10));
ylim([-1400 50]);
set(gca,'XTick',1:(model_num-1),'XTickLabel',model_short(2:end));
ylabel('log_{10} Bayes factor');
title('Bayes factor per subject');
subplot(2,1,2);
bar(rel_fit_bic(:,2:end)'/log(10));
ylim([-250 50]);
set(gca,'XTick',1:(model_num-1),'XTickLabel',model_short(2:end));
xlabel('BIC - BIC optim');
%legend(model_short);

% average llh fits (with SEM)
figure;  hold on;
errorbar(2:model_num,mean(rel_fit_llh(:,2:end),1),...
    sqrt(var(rel_fit_llh(:,2:end),[],1)/subj_num),'ko-','LineWidth',1.5);
plot(get(gca,'XLim'),[0 0], 'k--');
set(gca,'XTick',2:model_num,'XTickLabel',model_short(2:end));
ylabel('LLH - LLH optim');

% average bic fits (with SEM)
figure;
subplot(2,1,1); hold on;
bar(mean(rel_fit_bic(:,2:end)/log(10),1),'EdgeColor','none','FaceColor',[1 1 1]*0.5);
errorbar(1:(model_num-1),mean(rel_fit_bic(:,2:end)/log(10),1),...
    sqrt(var(rel_fit_bic(:,2:end)/log(10),[],1)/subj_num),'ko','LineWidth',1.5);
plot(get(gca,'XLim'),[-2 -2], 'k--');
ylim([-400 10]);
set(gca,'XTick',1:(model_num-1),'XTickLabel',model_short(2:end));
ylabel('log_{10} Bayes factor, mean +/- SEM');
title('avg. Bayes factor');
subplot(2,1,2); hold on;
bar(mean(rel_fit_bic(:,2:end)/log(10),1),'EdgeColor','none','FaceColor',[1 1 1]*0.5);
errorbar(1:(model_num-1),mean(rel_fit_bic(:,2:end)/log(10),1),...
    sqrt(var(rel_fit_bic(:,2:end)/log(10),[],1)/subj_num),'ko','LineWidth',1.5);
plot(get(gca,'XLim'),[-2 -2], 'k--');
ylim([-50 10]);

% normal and adjusted corefficient of determination
figure; hold on;
bar([fit_r2(:,1) fit_r2adj(:,1)]);
ylim([0.8 1]); xlim([0 (subj_num+1)]);
legend('R^2','R_{adj}^2', 'Location','SouthEast');
set(gca,'XTick',1:subj_num,'XTickLabel',subj_ids);
ylabel('R^2 and adjusted R^2');
title(sprintf('R2 and adjusted R2 for %s model', model_ids{1}),'Interpreter','None');

% sum of relative log-likelihoods
figure; hold on;
subplot(2,1,1);
bar(sum(rel_fit_llh(:,2:end),1));
ylabel(sprintf('LLH(model) - LLH(%s)', model_ids{1}),'Interpreter','None');
xlim([0 model_num]);
set(gca,'XTick',[],'XTickLabel',[]);
subplot(2,1,2);
bar(sum(rel_fit_bic(:,2:end),1)/log(10));
xlim([0 model_num]);
set(gca,'XTick',1:(model_num-1),'XTickLabel',model_short(2:end));
ylabel(sprintf('Model evidence (log10) model - %s', model_ids{1}),'Interpreter','None');

%figure;
%barh(sum(rel_fit_bic(:,2:end),1)/log(10));
%set(gca,'YTick',1:(model_num-1),'YTickLabel',model_short(2:end));


%% fancy plots for paper
% adjusted R2
pbar = 4/3;
figure('Color','white');
bar(fit_r2adj(:,1),'FaceColor',[1 1 1]*0.5,'EdgeColor','none');
xlim([0.5 (subj_num+0.5)]);  ylim([0 1]);
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
    'XTick',1:subj_num,'XTickLabel',{},'YTick',0:0.1:1);
xlabel('subject');
ylabel('adjusted R2');

% relative bic
pbar = 4/3;
figure('Color','white');
barh(-sum(rel_fit_bic(:,2:end),1)/log(10),'FaceColor',[1 1 1]*0.5,'EdgeColor','none');
xlim([0 3500]); ylim([0.5 (model_num-0.5)]);
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
    'YTick',1:(model_num-1),'YTickLabel',model_short(2:end),'YDir','reverse');
figure('Color','white');
barh(-sum(rel_fit_bic(:,2:end),1)/log(10),'FaceColor',[1 1 1]*0.5,'EdgeColor','none');
xlim([0 500]); ylim([0.5 (model_num-0.5)]);
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
    'YTick',1:(model_num-1),'YTickLabel',model_short(2:end),'YDir','reverse');

% relative bic per subject
pbar = 1;
bwidth_model = 0.8;
bprcnt_subj = 0.7;
bwidthf_subj = bwidth_model / (subj_num-1+bprcnt_subj);
bwidth_subj = bprcnt_subj*bwidthf_subj;
for i = 1:2
    figure('Color','white');
    if i == 1, xlim([-10 610]); else xlim([-10 20]); end
    ylim([0.5 (model_num-0.5)]);
    for model_idx = 1:(model_num-1)
        for subj_idx = 1:subj_num
            x = 0; w = -rel_fit_bic(subj_idx,model_idx+1)/log(10);
            if w < 0,  x = w;  w = abs(w);  end
            y = model_idx+bwidth_model/2-(subj_idx-1)*bwidthf_subj-bwidth_subj;
            h = bwidth_subj;
            rectangle('Position',[x y w h],...
                'EdgeColor','none','FaceColor',subj_colors(subj_idx,:));
        end  
    end
    set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
        'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
        'YTick',1:(model_num-1),'YTickLabel',model_short(2:end));
end
function plot_bic_comparison
%% plots the bic comparison for some model fits


%% settings
subj_ids = {'231139','231140','548203','548214','950091','950107','967716','eliana','kalpana','jason'};
%model_ids = {'mf_optim', 'mf_sepk', 'mf_sepk_nobias', 'mf_sepk_onetnd', 'mf_sepk_nobias_onetnd'};
model_ids = {'mf_optim', 'mf_sepw', 'mf_sepw_nobias', 'mf_sepw_onetnd', 'mf_sepw_nobias_onetnd'};
subj_num = length(subj_ids);
model_num = length(model_ids);
data_folder = 'fit_data';
comb_model_idx = 2;


%% collect data
bics = zeros(subj_num, model_num);
llhs = zeros(subj_num, model_num);
for subj_idx = 1:subj_num
    for model_idx = 1:model_num
        % load data and process
        d = load([data_folder filesep ...
            model_ids{model_idx} '_' subj_ids{subj_idx} '.mat']);
        n = d.dstats.n;
        k = length(d.p_names);
        llhs(subj_idx,model_idx) = d.best_llh;
        bics(subj_idx,model_idx) = -d.best_llh + 0.5*k*log(n);
    end
end
rel_llhs = bsxfun(@minus,llhs,llhs(:,comb_model_idx));
rel_bics = bsxfun(@minus,bics,bics(:,comb_model_idx));



%% plot full bic and likelihood
figure;
subplot(2,1,1);
bar(llhs);
set(gca,'Layer','top','TickDir','out','FontName','Helvetica','FontSize',12,...
    'XTick',1:subj_num,'XTickLabel',{});
ylabel('LLH');
legend(model_ids);
subplot(2,1,2);
bar(bics);
set(gca,'Layer','top','TickDir','out','FontName','Helvetica','FontSize',12,...
    'XTick',1:subj_num,'XTickLabel',subj_ids);
ylabel('BIC');


%% plot relativ bic and likelihood
figure;
subplot(2,1,1);
bar(rel_llhs);
set(gca,'Layer','top','TickDir','out','FontName','Helvetica','FontSize',12,...
    'XTick',1:subj_num,'XTickLabel',{});
ylabel(['LLH - LLH(' model_ids{comb_model_idx} ')']);
legend(model_ids);
subplot(2,1,2);
bar(rel_bics);
set(gca,'Layer','top','TickDir','out','FontName','Helvetica','FontSize',12,...
    'XTick',1:subj_num,'XTickLabel',subj_ids);
ylabel(['BIC - BIC(' model_ids{comb_model_idx} ')']);


%% average comparison
figure;
subplot(2,1,1);  hold on;
bar(mean(rel_llhs,1),1);
for model_idx = 1:model_num
    s = [mean(rel_llhs(:,model_idx)) sqrt(var(rel_llhs(:,model_idx))/subj_num)];
    if s(2) > 0
        plot(ones(1,2)*model_idx, [(s(1)-s(2)) (s(1)+s(2))], 'k-');
    end
end
set(gca,'Layer','top','TickDir','out','FontName','Helvetica','FontSize',12,...
    'XTick',1:model_num,'XTickLabel',{});
ylabel('avg LLH');

subplot(2,1,2);  hold on;
bar(mean(rel_bics,1),1);
for model_idx = 1:model_num
    s = [mean(rel_bics(:,model_idx)) sqrt(var(rel_bics(:,model_idx))/subj_num)];
    if s(2) > 0
        plot(ones(1,2)*model_idx, [(s(1)-s(2)) (s(1)+s(2))], 'k-');
    end
end
set(gca,'Layer','top','TickDir','out','FontName','Helvetica','FontSize',12,...
    'XTick',1:model_num,'XTickLabel',model_ids);
ylabel('avg BIC');

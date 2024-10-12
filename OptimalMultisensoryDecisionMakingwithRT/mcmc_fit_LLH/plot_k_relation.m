function plot_k_relation
%% plots the relation between the fitted sensitivities


fit_folder = 'fit_data';
var_type = 0;       % 0 - none, 1 - powsep, 2 - powcomb
burnin = 0;% 22000;
%burninw = 0;
conf_lo = 2.5;
conf_up = 97.5;
pred_samples = 1000;
coh_range = [-0.1 0.8];
cohs = [0.001 0.12 0.25 0.37 0.52 0.7];
%cohs = [0.25 0.37 0.7];
coh_num = length(cohs);

sepk_fit_file = @(subj_id) [fit_folder filesep 'mf_sepk_' subj_id '.mat'];
sepw_fit_file = @(subj_id) [fit_folder filesep 'mf_sepw_' subj_id '.mat'];
sepk_fit_file_powsep = @(subj_id) [fit_folder filesep 'mf_sepk_powsep_' subj_id '.mat'];
sepw_fit_file_powsep = @(subj_id) [fit_folder filesep 'mf_sepw_powsep_' subj_id '.mat'];
sepk_fit_file_powcomb = @(subj_id) [fit_folder filesep 'mf_sepk_powcomb_' subj_id '.mat'];
sepw_fit_file_powcomb = @(subj_id) [fit_folder filesep 'mf_sepw_powcomb_' subj_id '.mat'];

%subj_ids = {'231139','231140','548203','548214','950091','950107','967716','eliana','kalpana','jason','kalpana_cut'};
subj_ids = {'231139','231140','548203','548214','950091','950107','967716','jason','eliana','kalpana'};
%subj_ids = {'231139','231140','548203','548214','950091','950107','967716'};
%subj_ids = {'231139','231140','548203','548214','950091','950107','967716','eliana','kalpana_cut','jason'};


% compute unimodal/combined k's (with stats) from samples
subj_num = length(subj_ids);
sstats = cell(1, subj_num);
ravg_est_w = NaN(subj_num,coh_num,7);
ravg_sstats = zeros(subj_num,coh_num,4,7);
avg_sstats = cell(coh_num,7);
for subj_idx = 1:subj_num
    % load file & pre-process file
    if var_type == 0
        % sepk model - everything separate
        bound_name = 'free bound';
        d = load(sepk_fit_file(subj_ids{subj_idx}));
        dw = load(sepw_fit_file(subj_ids{subj_idx}));
        cs = d.cohs;
        c_num = length(cs);
        s = d.s((burnin+1):end,:);
        sw = dw.s((burnin+1):end,:);
        kvis_idx = 0;
        kvest_idx = 2*c_num+1;
        kcomb_idx = 2*c_num+2;
        best_p = d.best_p;
        best_pw = dw.best_p;
    elseif var_type == 1
        % sepk_powsep model, with parameters
        % 1-avis, 2-bvis, 3-gvis, 4-gam1vis, 5-gam2vis, 6-kvest, 7-boundvest,
        % 8-kcomb1, 9-kcomb2, ...
        % using
        % k_vis(c) = a_vis c^gam1_vis / sqrt(1+b_vis c^gam2_vis)
        bound_name = 'powsep bound';
        d = load(sepk_fit_file_powsep(subj_ids{subj_idx}));
        %dw = load(sepw_fit_file_powsep(subj_ids{subj_ids}));
        cs = d.cohs;
        c_num = length(cs);
        s = d.s((burnin+1):end,:);
        % compute kvis samples and assemble new sample matrix
        kvis_s = bsxfun(@times,s(:,1),bsxfun(@power,cs,s(:,4)))./...
            sqrt(1+bsxfun(@times,s(:,2),bsxfun(@power,cs,s(:,5))));
        s = [kvis_s s(:,6) s(:,8:(8+c_num-1))];
        kvis_idx = 0;
        kvest_idx = c_num+1;
        kcomb_idx = c_num+1;
        % do the same with best_p
        best_p = d.best_p;
        best_p = [(best_p(1)*cs.^best_p(4)./sqrt(1+best_p(2)*cs.^best_p(5))) ...
            best_p(6) best_p(8:(8+c_num-1))];
    elseif var_type == 2
        % sepk_powcomb model, with parameters
        % 1-avis, 2-bvis, 3-gvis, 4-gamvis, 5-kvest, 6-boundvest,
        % 7-kcomb1, 8-kcomb2, ...
        % using
        % k_vis(c) = a_vis c^gam_vis / sqrt(1+b_vis c^gam_vis)
        bound_name = 'powcomb bound';
        d = load(sepk_fit_file_powcomb(subj_ids{subj_idx}));
        %dw = load(sepw_fit_file_powcomb(subj_ids{subj_idx}));
        cs = d.cohs;
        c_num = length(cs);
        s = d.s((burnin+1):end,:);
        % compute kvis samples and assemble new sample matrix
        kvis_s = bsxfun(@times,s(:,1),bsxfun(@power,cs,s(:,4)))./...
            sqrt(1+bsxfun(@times,s(:,2),bsxfun(@power,cs,s(:,4))));
        s = [kvis_s s(:,5) s(:,7:(7+c_num-1))];
        kvis_idx = 0;
        kvest_idx = c_num+1;
        kcomb_idx = c_num+1;
        % do the same with best_p
        best_p = d.best_p;
        best_p = [(best_p(1)*cs.^best_p(4)./sqrt(1+best_p(2)*cs.^best_p(4))) ...
            best_p(5) best_p(7:(7+c_num-1))];
    else
        error('Unknown var_type');
    end
    s_num = size(s,1);
    % compute parameter mean / confidence intervals
    kvest = [mean(s(:,kvest_idx)) ...
             prctile(s(:,kvest_idx),conf_lo) ...
             prctile(s(:,kvest_idx),conf_up) ...
             best_p(kvest_idx)];
    kvis = zeros(c_num,4);   % mean lo up mode
    kcomb = zeros(c_num,4);
    kpred = zeros(c_num,4);
    kpred_nw = zeros(c_num,4);
    wvis = zeros(c_num,4);
    wvis_pred = zeros(c_num,4);
    for c_idx = 1:c_num
        % compute predicted kcomb by subsampling (with replacement)
        kvis_s2 = s(floor(s_num*rand(pred_samples,1))+1,kvis_idx+c_idx).^2;
        kvest_s2 = s(floor(s_num*rand(pred_samples,1))+1,kvest_idx).^2;
        kpred_s = sqrt(kvis_s2+kvest_s2);
        kpred_nw_s = sqrt(4 * kvis_s2 .* kvest_s2 ./ (kvis_s2 + kvest_s2));
        % percentiles of these 
        kvis(c_idx,:) = [mean(s(:,kvis_idx+c_idx)) ...
            prctile(s(:,kvis_idx+c_idx),conf_lo) ...
            prctile(s(:,kvis_idx+c_idx),conf_up) ...
            best_p(kvis_idx+c_idx)];
        kcomb(c_idx,:) = [mean(s(:,kcomb_idx+c_idx)) ...
            prctile(s(:,kcomb_idx+c_idx),conf_lo) ...
            prctile(s(:,kcomb_idx+c_idx),conf_up) ...
            best_p(kcomb_idx+c_idx)];
        kpred(c_idx,:) = [mean(kpred_s) ...
            prctile(kpred_s,conf_lo) prctile(kpred_s,conf_up)...
            sqrt(best_p(kvest_idx)^2+best_p(kvis_idx+c_idx)^2)];
        kpred_nw(c_idx,:) = [mean(kpred_nw_s) ...
            prctile(kpred_nw_s,conf_lo) prctile(kpred_nw_s,conf_up)...
            sqrt(4 * best_p(kvest_idx)^2 * best_p(kvis_idx+c_idx)^2 / ...
                 (best_p(kvest_idx)^2+best_p(kvis_idx+c_idx)^2))];
        % combination weight statistics
        %w = s(:,kvis_idx+c_idx).^2./(s(:,kvis_idx+c_idx).^2+s(:,kvest_idx));
        %wvis(c_idx,:) = [mean(w) ...
        %    prctile(w,conf_lo) prctile(w,conf_up)...
        %    (kvis(c_idx,4)^2/(kvis(c_idx,4)^2+kvest(4)^2))];
        w = sw(:,kvis_idx+c_idx).^2./(sw(:,kvis_idx+c_idx).^2+sw(:,kvest_idx).^2);
        wvis(c_idx,:) = [mean(sw(:,kcomb_idx+c_idx)) ...
            prctile(sw(:,kcomb_idx+c_idx),conf_lo) ...
            prctile(sw(:,kcomb_idx+c_idx),conf_up) ...
            best_pw(kcomb_idx+c_idx)];
        wvis_pred(c_idx,:) = [mean(w) ...
            prctile(w,conf_lo) prctile(w,conf_up)...
            (kvis(c_idx,4)^2/(kvis(c_idx,4)^2+kvest(4)^2))];
        % compute estimator variance for robust averaging
        all_c_idx = find(abs(cs(c_idx) - cohs) < 1e-5);
        if isempty(all_c_idx), error('Unknown coherence value'); end
        % using Fisher transform to get variance of visual weights
        ravg_est_w(subj_idx,all_c_idx,:) = ...
            [var(s(:,kcomb_idx+c_idx)) var(kpred_s) var(kpred_nw_s) ...
             var(atanh(2*sw(:,kcomb_idx+c_idx)-1)) var(atanh(2*w-1)) ...
             var(s(:,kvis_idx+c_idx)) var(s(:,kvest_idx))];
        ravg_sstats(subj_idx,all_c_idx,:,1) = kcomb(c_idx,:);
        ravg_sstats(subj_idx,all_c_idx,:,2) = kpred(c_idx,:);
        ravg_sstats(subj_idx,all_c_idx,:,3) = kpred_nw(c_idx,:);
        ravg_sstats(subj_idx,all_c_idx,:,4) = wvis(c_idx,:);
        ravg_sstats(subj_idx,all_c_idx,:,5) = wvis_pred(c_idx,:);
        ravg_sstats(subj_idx,all_c_idx,:,6) = kvis(c_idx,:);
        ravg_sstats(subj_idx,all_c_idx,:,7) = kvest;
        avg_sstats{all_c_idx,1} = cat(2,avg_sstats{all_c_idx,1},kcomb(c_idx,4));
        avg_sstats{all_c_idx,2} = cat(2,avg_sstats{all_c_idx,2},kpred(c_idx,4));
        avg_sstats{all_c_idx,3} = cat(2,avg_sstats{all_c_idx,3},kpred_nw(c_idx,4));
        avg_sstats{all_c_idx,4} = cat(2,avg_sstats{all_c_idx,4},wvis(c_idx,4));
        avg_sstats{all_c_idx,5} = cat(2,avg_sstats{all_c_idx,5},wvis_pred(c_idx,4));
        avg_sstats{all_c_idx,6} = cat(2,avg_sstats{all_c_idx,6},kvis(c_idx,4));
        avg_sstats{all_c_idx,7} = cat(2,avg_sstats{all_c_idx,7},kvest(4));
    end
    % create data structure to save
    sstats{subj_idx} = struct('cohs',cs,...
        'kvis',kvis,'kvest',kvest','kcomb',kcomb,'kpred',kpred,'kpred_nw',kpred_nw, ...
        'wvis',wvis,'wvis_pred',wvis_pred);
end

% average over subject by normal averaging + SEM
avg_kcomb = zeros(2,coh_num);
avg_kpred = zeros(2,coh_num);
avg_kpred_nw = zeros(2,coh_num);
avg_wvis = zeros(2,coh_num);
avg_wvis_pred = zeros(2,coh_num);
avg_diff = zeros(1,coh_num);
avg_diff_nw = zeros(1,coh_num);
avg_diff_wvis = zeros(1,coh_num);
avg_kvis = zeros(2,coh_num);
avg_kvest = zeros(2,coh_num);
for coh_idx = 1:coh_num
    avg_kcomb(:,coh_idx) = [mean(avg_sstats{coh_idx,1}) ...
        sqrt(var(avg_sstats{coh_idx,1})/length(avg_sstats{coh_idx,1}))];
    avg_kpred(:,coh_idx) = [mean(avg_sstats{coh_idx,2}) ...
        sqrt(var(avg_sstats{coh_idx,2})/length(avg_sstats{coh_idx,2}))];
    avg_kpred_nw(:,coh_idx) = [mean(avg_sstats{coh_idx,3}) ...
        sqrt(var(avg_sstats{coh_idx,3})/length(avg_sstats{coh_idx,3}))];
    avg_wvis(:,coh_idx) = [mean(avg_sstats{coh_idx,4}) ...
        sqrt(var(avg_sstats{coh_idx,4})/length(avg_sstats{coh_idx,4}))];
    avg_wvis_pred(:,coh_idx) = [mean(avg_sstats{coh_idx,5}) ...
        sqrt(var(avg_sstats{coh_idx,5})/length(avg_sstats{coh_idx,5}))];
    avg_kvis(:,coh_idx) = [mean(avg_sstats{coh_idx,6}) ...
        sqrt(var(avg_sstats{coh_idx,6})/length(avg_sstats{coh_idx,6}))];
    avg_kvest(:,coh_idx) = [mean(avg_sstats{coh_idx,7}) ...
        sqrt(var(avg_sstats{coh_idx,7})/length(avg_sstats{coh_idx,7}))];
    [~,avg_diff(coh_idx)] = ttest(avg_sstats{coh_idx,1},avg_sstats{coh_idx,2});
    [~,avg_diff_nw(coh_idx)] = ttest(avg_sstats{coh_idx,1},avg_sstats{coh_idx,3});
    [~,avg_diff_wvis(coh_idx)] = ttest(avg_sstats{coh_idx,4},avg_sstats{coh_idx,5});
end


% average over subjects by robust averaging
no_data_idx = isnan(ravg_est_w);
ravg_est_w = 1 ./ ravg_est_w;
ravg_est_w = bsxfun(@rdivide,ravg_est_w,nansum(ravg_est_w,1));
ravg_est_w(no_data_idx) = 0;
ravg_kcomb = squeeze(sum(bsxfun(@times,ravg_est_w(:,:,1),...
                                       ravg_sstats(:,:,:,1)),1)); 
ravg_kpred = squeeze(sum(bsxfun(@times,ravg_est_w(:,:,2),...
                                       ravg_sstats(:,:,:,2)),1));
ravg_kpred_nw = squeeze(sum(bsxfun(@times,ravg_est_w(:,:,3),...
                                       ravg_sstats(:,:,:,3)),1));
ravg_wvis = squeeze(sum(bsxfun(@times,ravg_est_w(:,:,4),...
                                       ravg_sstats(:,:,:,4)),1));                                 
ravg_wvis_pred = squeeze(sum(bsxfun(@times,ravg_est_w(:,:,5),...
                                           ravg_sstats(:,:,:,5)),1));                                 

                                       
% plot results
for subj_idx = 1:subj_num
    figure;
    subplot(2,1,1); hold on;
    s = sstats{subj_idx};
    errorbar(s.cohs-0.01,s.kpred_nw(:,4),...
        s.kpred_nw(:,4)-s.kpred_nw(:,2),s.kpred_nw(:,3)-s.kpred_nw(:,4),...
        'd-','Color',[0.8 0 0],'MarkerSize',10,'MarkerFaceColor',[1 1 1]);        
    errorbar(s.cohs+0.01,s.kpred(:,4),...
        s.kpred(:,4)-s.kpred(:,2),s.kpred(:,3)-s.kpred(:,4),...
        'o-','Color',[0 0.8 0],'MarkerSize',10,'MarkerFaceColor',[1 1 1]);
    errorbar(s.cohs,s.kcomb(:,4),...
        s.kcomb(:,4)-s.kcomb(:,2),s.kcomb(:,3)-s.kcomb(:,4),...
        'ks','MarkerSize',10,'MarkerFaceColor',[1 1 1]);
    % collect significant differences
    sig_diff_cs = s.cohs(...
        s.kpred_nw(:,2)>s.kcomb(:,3) | s.kcomb(:,2)>s.kpred_nw(:,3));
    if ~isempty(sig_diff_cs)
        ylims = get(gca,'YLim');
        plot(sig_diff_cs+0.01,...
            ones(size(sig_diff_cs))*(ylims(1)+0.05*(ylims(2)-ylims(1))),...
            '*','Color',[0.8 0 0],'MarkerSize',10);
    end
    sig_diff_cs = s.cohs(...
        s.kpred(:,2)>s.kcomb(:,3) | s.kcomb(:,2)>s.kpred(:,3));
    if ~isempty(sig_diff_cs)
        ylims = get(gca,'YLim');
        plot(sig_diff_cs+0.01,...
            ones(size(sig_diff_cs))*(ylims(1)+0.05*(ylims(2)-ylims(1))),...
            '*','Color',[0 0.8 0],'MarkerSize',10);
    end
    % label graphs
    xlim(coh_range);
    ylabel('sensitivity');
    legend({'even weighting','corr. weighting','observed'},'Location','NorthWest');
    title(['subj ' subj_ids{subj_idx} ', ' bound_name]);
    
    subplot(2,1,2); hold on;
    
    for i = 1:length(s.cohs)
        h = plot(ones(1,2)*(s.cohs(i)+0.005),s.wvis_pred(i,2:3),...
            'k-','LineWidth',0.5);
        set(get(get(h,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off')
        h = plot(ones(1,2)*(s.cohs(i)-0.005),s.wvis(i,2:3),...
            'k-','LineWidth',0.5);
        set(get(get(h,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off')
    end
    plot(s.cohs+0.005,s.wvis_pred(:,1),...
        'ko','MarkerSize',10,'MarkerFaceColor',[1 1 1]);
    %errorbar(s.cohs+0.005,s.wvis_pred(:,4),...
    %    s.wvis_pred(:,4)-s.wvis_pred(:,2),s.wvis_pred(:,3)-s.wvis_pred(:,4),...
    %    'ko-','MarkerSize',10,'MarkerFaceColor',[1 1 1]);
    %errorbar(s.cohs-0.005,s.wvis(:,1),...
    %    s.wvis(:,1)-s.wvis(:,2),s.wvis(:,3)-s.wvis(:,1),...
    %    'ks','MarkerSize',10,'MarkerFaceColor',[1 1 1]);
    plot(s.cohs-0.005,s.wvis(:,1),'ks','MarkerSize',10,'MarkerFaceColor',[1 1 1]);
    % collect significant differences
    sig_diff_cs = s.cohs(...
        min(s.wvis_pred(:,2),s.wvis_pred(:,4))>max(s.wvis(:,3),s.wvis(:,4)) | ...
        min(s.wvis(:,2),s.wvis(:,4))>max(s.wvis_pred(:,3),s.wvis_pred(:,4)));
    if ~isempty(sig_diff_cs)
        ylims = get(gca,'YLim');
        plot(sig_diff_cs,...
            ones(size(sig_diff_cs))*(ylims(1)+0.05*(ylims(2)-ylims(1))),...
            'k*','MarkerSize',10);
    end
    xlim(coh_range);  ylim([0 1]);
    ylabel('vis weight');
    xlabel('coherence');
    legend({'pred','observed'},'Location','NorthWest');
end


% fancy-plot average
pbar = 7/6;
figure('Color','white');  hold on;
xlim([-0.02 0.72]);
sd_mul = 1.96;
errorshade(cohs,avg_kvis(1,:)-sd_mul*avg_kvis(2,:),avg_kvis(1,:)+sd_mul*avg_kvis(2,:),0.5*[0 0 200]/255+0.5*[1 1 1]);
errorshade(cohs,avg_kvest(1,:)-sd_mul*avg_kvest(2,:),avg_kvest(1,:)+sd_mul*avg_kvest(2,:),0.5*[0 200 0]/255+0.5*[1 1 1]);
errorshade(cohs,avg_kpred(1,:)-sd_mul*avg_kpred(2,:),avg_kpred(1,:)+sd_mul*avg_kpred(2,:),[200 0 0]/255);
plot(cohs,avg_kvis(1,:),'LineWidth',1.5,'Color',0.5*[0 0 200]/255+0.5*[1 1 1]);
plot(cohs,avg_kvest(1,:),'LineWidth',1.5,'Color',0.5*[0 200 0]/255+0.5*[1 1 1]);
plot(cohs,avg_kvis(1,:),'o','Color',0.5*[0 0 200]/255+0.5*[1 1 1],'MarkerSize',4,'MarkerFaceColor',[1 1 1]);
plot(cohs,avg_kvest(1,:),'o','Color',0.5*[0 200 0]/255+0.5*[1 1 1],'MarkerSize',4,'MarkerFaceColor',[1 1 1]);
plot(cohs,avg_kpred(1,:),'LineWidth',1.5,'Color',[200 0 0]/255);
plot(cohs,avg_kpred(1,:),'o','Color',[200 0 0]/255,'MarkerSize',4,'MarkerFaceColor',[1 1 1]);
for c_idx = 1:length(cohs)
    plot(cohs(c_idx)*[1 1]+0.01,avg_kcomb(1,c_idx)+sd_mul*avg_kcomb(2,c_idx)*[-1 1],'Color',[200 0 0]/255);
end
plot(cohs+0.01,avg_kcomb(1,:),'s','Color',[200 0 0]/255,'MarkerSize',4,'MarkerFaceColor',[1 1 1]);
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
    'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
    'XTick',cohs,'XTickLabel',{'0','12','25','37','51','70'});
xlabel('coherence');
ylabel('sensitivity');


% plot average
figure;
%subplot(2,1,1);
hold on;
%errorshade(cohs-0.01,avg_kpred_nw(1,:)-avg_kpred_nw(2,:),...
%    avg_kpred_nw(1,:)+avg_kpred_nw(2,:),[0.8 0 0]);
errorshade(cohs+0.01,avg_kpred(1,:)-avg_kpred(2,:),...
    avg_kpred(1,:)+avg_kpred(2,:),[0 0.8 0]);
%plot(cohs-0.01,avg_kpred_nw(1,:),...
%    'd-','Color',[0.8 0 0],'MarkerSize',10,'MarkerFaceColor',[1 1 1],'LineWidth',1);
plot(cohs+0.01,avg_kpred(1,:),...
    'o-','Color',[0 0.8 0],'MarkerSize',10','MarkerFaceColor',[1 1 1],'LineWidth',1);
errorbar(cohs,avg_kcomb(1,:),avg_kcomb(2,:),...
    'ks','MarkerSize',10,'MarkerFaceColor',[1 1 1],'LineWidth',1);

% mark significant differences
%sig_diff_cs = cohs(avg_diff_nw < 0.05);
%if ~isempty(sig_diff_cs)
%    ylims = get(gca,'YLim');
%        plot(sig_diff_cs-0.01,...
%            ones(size(sig_diff_cs))*(ylims(1)+0.05*(ylims(2)-ylims(1))),...
%            '*','Color',[0.8 0 0],'MarkerSize',10);
%end
sig_diff_cs = cohs(avg_diff < 0.05);
if ~isempty(sig_diff_cs)
    ylims = get(gca,'YLim');
        plot(sig_diff_cs+0.01,...
            ones(size(sig_diff_cs))*(ylims(1)+0.05*(ylims(2)-ylims(1))),...
            '*','Color',[0 0.8 0],'MarkerSize',10);
end
xlim(coh_range);
ylabel('sensitivity');
%legend({'even weighting','corr. weighting','observed'},'Location','NorthWest');
legend({'corr. weighting','observed'},'Location','NorthWest');
title(['all subjects, average of mode + SEM, ' bound_name]);
figure;
%subplot(2,1,2);
hold on;
errorshade(cohs+0.005,avg_wvis_pred(1,:)-avg_wvis_pred(2,:),...
    avg_wvis_pred(1,:)+avg_wvis_pred(2,:));
plot(cohs+0.005,avg_wvis_pred(1,:),...
    'ko-','MarkerSize',10','MarkerFaceColor',[1 1 1],'LineWidth',1);
errorbar(cohs-0.005,avg_wvis(1,:),avg_wvis(2,:),...
    'ks','MarkerSize',10,'MarkerFaceColor',[1 1 1],'LineWidth',1);
% mark significant differences
sig_diff_cs = cohs(avg_diff_wvis < 0.05);
if ~isempty(sig_diff_cs)
    ylims = get(gca,'YLim');
        plot(sig_diff_cs,...
            ones(size(sig_diff_cs))*(ylims(1)+0.05*(ylims(2)-ylims(1))),...
            'k*','MarkerSize',10);
end
xlim(coh_range);  ylim([0 1]);
ylabel('vis weight');
xlabel('coherence');
legend({'pred','observed'},'Location','NorthWest');


% plot robust average
figure;
subplot(2,1,1);  hold on;
errorshade(cohs-0.01,ravg_kpred_nw(:,2),ravg_kpred_nw(:,3),[0.8 0 0]);
errorshade(cohs+0.01,ravg_kpred(:,2),ravg_kpred(:,3),[0 0.8 0]);
plot(cohs-0.01,ravg_kpred_nw(:,4),...
    'd-','Color',[0.8 0 0],'MarkerSize',10,'MarkerFaceColor',[1 1 1],'LineWidth',1);
plot(cohs+0.01,ravg_kpred(:,4),...
    'o-','Color',[0 0.8 0],'MarkerSize',10','MarkerFaceColor',[1 1 1],'LineWidth',1);
errorbar(cohs,ravg_kcomb(:,4),...
    ravg_kcomb(:,4)-ravg_kcomb(:,2),ravg_kcomb(:,3)-ravg_kcomb(:,4),...
    'ks','MarkerSize',10','MarkerFaceColor',[1 1 1],'LineWidth',1);
% mark significant differences
sig_diff_cs = cohs(...
     ravg_kpred_nw(:,2)>ravg_kcomb(:,3) | ravg_kcomb(:,2)>ravg_kpred_nw(:,3));
if ~isempty(sig_diff_cs)
    ylims = get(gca,'YLim');
    plot(sig_diff_cs-0.01,...
        ones(size(sig_diff_cs))*(ylims(1)+0.05*(ylims(2)-ylims(1))),...
        '*','Color',[0.8 0 0],'MarkerSize',10);
end
sig_diff_cs = cohs(...
     ravg_kpred(:,2)>ravg_kcomb(:,3) | ravg_kcomb(:,2)>ravg_kpred(:,3));
if ~isempty(sig_diff_cs)
    ylims = get(gca,'YLim');
    plot(sig_diff_cs+0.01,...
        ones(size(sig_diff_cs))*(ylims(1)+0.05*(ylims(2)-ylims(1))),...
        '*','Color',[0 0.8 0],'MarkerSize',10);
end
xlim(coh_range);
ylabel('sensitivity');
legend({'even weighting','corr. weighting','observed'},'Location','NorthWest');
title(['all subjects, robust averaging, ' bound_name]);
subplot(2,1,2);  hold on;
errorshade(cohs+0.005,ravg_wvis_pred(:,2),ravg_wvis_pred(:,3));
plot(cohs+0.005,ravg_wvis_pred(:,1),...
    'ko-','MarkerSize',10','MarkerFaceColor',[1 1 1],'LineWidth',1);
errorbar(cohs-0.005,ravg_wvis(:,1),ravg_wvis(:,1)-ravg_wvis(:,2),...
    ravg_wvis(:,3)-ravg_wvis(:,1),...
    'ks','MarkerSize',10','MarkerFaceColor',[1 1 1],'LineWidth',1);
xlim(coh_range);  ylim([0 1]);
xlabel('coherence');
ylabel('vis weight');
legend({'pred','observed'},'Location','NorthWest');




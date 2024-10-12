function plot_fits(subj_id, model_fn)
%% plot currently best fits for given subject and model

%% settings
vis_base_cols = [200 0 0; 200 180 0] / 255;
vest_base_col = [0 0 200] / 255;
comb_base_cols = [0 200 0; 0 180 200] / 255;    % changed green and blue according to my habit
pbar = 8/3;
p_conf = 0.95;
rt_lims = [0.6 1.4];
plot_by_coh = true;
plot_by_mod = false;
plot_incorr_rts = false;
% data_folder = 'fit_data';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======= Change the file directory===========

data_folder = 'D:\Paper_rawdata\Code\Labtools\OptimalMultisensoryDecisionMakingwithRT-main\data\vest_vis_gauss_ZZ\fit_data';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isstr(model_fn), model_name = model_fn;
else model_name = func2str(model_fn); end
data_filename = [data_folder filesep model_name '_' subj_id '.mat'];


%% load model fits, and set up colors
if ~exist(data_filename,'file')
    error('Fit data file "%s" not found', data_filename);
end
d = load(data_filename);
if isempty(d.mstats.vis), vis_num = 0;
else
    vis_num = length(d.mstats.vis);
    % special case of 1, 3 or 6 visual conditions
    if vis_num == 1, vis_cols = vis_base_cols(1,:);
    elseif vis_num == 3
        vis_cols = [0  50 200; 0 100 200; 0 180 200] / 255;
    elseif vis_num == 6
        vis_cols = [0   0 200; 0  25 200; 0  50 200; ...
                    0 100 200; 0 140 200; 0 180 200] / 255;
    else
        vis_cols = bsxfun(@plus, vis_base_cols(1,:), ...
            bsxfun(@times, (vis_base_cols(2,:)-vis_base_cols(1,:)), ...
                (0:(vis_num-1))'/(vis_num-1)));
    end
end
if isempty(d.mstats.vest), vest_num = 0;
else
    vest_num = 1;
    vest_col = vest_base_col;
end
if isempty(d.mstats.comb), comb_num = 0;
else
    comb_num = length(d.mstats.comb);
    if comb_num == 1, comb_cols = comb_base_cols(1,:);
    elseif comb_num == 3
        comb_cols = [200  50 0; 200 100 0; 200 180 0] / 255;
    elseif comb_num == 6
        comb_cols = [200   0 0; 200  25 0; 200  50 0; ...
                     200 100 0; 200 140 0; 200 180 0] / 255;
    else
        comb_cols = bsxfun(@plus, comb_base_cols(1,:), ...
            bsxfun(@times, (comb_base_cols(2,:)-comb_base_cols(1,:)), ...
                (0:(comb_num-1))'/(comb_num-1)));
    end
end


%% load subject data to get error bars on p(choice)
if ~isempty(p_conf)
    s = load_condi_data(subj_id,false);
    for c_idx = 1:vis_num
        sc = bin_by_es(s.vis(c_idx).choice,s.vis(c_idx).rt,s.vis(c_idx).es,...
                       'pconf',p_conf);
        d.dstats.vis(c_idx).pr_l = sc.pr_l;
        d.dstats.vis(c_idx).pr_u = sc.pr_u;
    end
    if vest_num > 0
        sc = bin_by_es(s.vest(1).choice,s.vest(1).rt,s.vest(1).es,'pconf',p_conf);
        d.dstats.vest(1).pr_l = sc.pr_l;
        d.dstats.vest(1).pr_u = sc.pr_u;
    end
    for c_idx = 1:comb_num
        sc = bin_by_es(s.comb(c_idx).choice,s.comb(c_idx).rt,s.comb(c_idx).es,...
                       'pconf',p_conf);
        d.dstats.comb(c_idx).pr_l = sc.pr_l;
        d.dstats.comb(c_idx).pr_u = sc.pr_u;
    end
end


%% plot by coherence
if plot_by_coh
    assert(vis_num == comb_num);
    for c_idx = 1:vis_num
        figure('Color','white');
        % chornometric
        subplot(2,1,1);  hold on;
        xlim([-17 17]);  if ~isempty(rt_lims), ylim(rt_lims); end
        if plot_incorr_rts
            % model
            plot(d.dstats.vis(c_idx).hs, d.mstats.vis(c_idx).rt_incorr,...
                '--','LineWidth',0.5,'Color',vis_cols(c_idx,:));
            if vest_num > 0
                plot(d.dstats.vest.hs, d.mstats.vest.rt_incorr,...
                    '--','LineWidth',0.5,'Color',vest_col);
            end
            plot(d.dstats.comb(c_idx).hs, d.mstats.comb(c_idx).rt_incorr,...
                '--','LineWidth',0.5,'Color',comb_cols(c_idx,:));
            % data
            plot_rt_incorr(d.dstats.vis(c_idx),vis_cols(c_idx,:));
            if vest_num > 0, plot_rt_incorr(d.dstats.vest,vest_col); end
            plot_rt_incorr(d.dstats.comb(c_idx),comb_cols(c_idx,:));
        end
        % model
        plot(d.dstats.vis(c_idx).hs, d.mstats.vis(c_idx).rt_corr, ...
            '-','LineWidth',1.5,'Color',vis_cols(c_idx,:));
        plot(d.dstats.vest.hs, d.mstats.vest.rt_corr, ...
            '-','LineWidth',1.5,'Color',vest_col);
        plot(d.dstats.comb(c_idx).hs, d.mstats.comb(c_idx).rt_corr, ...
            '-','LineWidth',1.5,'Color',comb_cols(c_idx,:));
        % data
        plot_rt_corr(d.dstats.vis(c_idx),vis_cols(c_idx,:));
        if vest_num > 0, plot_rt_corr(d.dstats.vest,vest_col); end
        plot_rt_corr(d.dstats.comb(c_idx),comb_cols(c_idx,:));
        % formatting
        set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
            'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
            'XTick',[-10,0,10],'XTickLabel',{});
        if ~isempty(rt_lims), set(gca,'YTick',rt_lims(1):0.2:rt_lims(2)); end
        ylabel('reaction time');
        
        % psychometric
        subplot(2,1,2);  hold on;
        xlim([-17 17]);  ylim([0 1]);
        % model
        plot(d.dstats.vis(c_idx).hs, d.mstats.vis(c_idx).pr, ...
            '-','LineWidth',1.5,'Color',vis_cols(c_idx,:));
        plot(d.dstats.vest.hs, d.mstats.vest.pr, ...
            '-','LineWidth',1.5,'Color',vest_col);
        plot(d.dstats.comb(c_idx).hs, d.mstats.comb(c_idx).pr, ...
            '-','LineWidth',1.5,'Color',comb_cols(c_idx,:));
        % data
        plot_pr(d.dstats.vis(c_idx),vis_cols(c_idx,:));
        if vest_num > 0, plot_pr(d.dstats.vest,vest_col); end
        plot_pr(d.dstats.comb(c_idx),comb_cols(c_idx,:));
        % formatting
        set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
            'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
            'XTick',[-10,0,10],'XTickLabel',{},...
            'YTick',0:025:1,'YTickLabel',{'0','','0.5','','1'});
        xlabel('heading');  ylabel('p(right)');

    end
end


%% plot by modality
if plot_by_mod
    % visual
    if vis_num > 0
        figure('Color','white');
        subplot(2,1,1);  hold on;
        xlim([-17 17]);
        if plot_incorr_rts
            for c_idx = 1:vis_num
                plot(d.dstats.vis(c_idx).hs, d.mstats.vis(c_idx).rt_incorr,...
                    '--','LineWidth',0.5,'Color',vis_cols(c_idx,:));
                plot_rt_incorr(d.dstats.vis(c_idx),vis_cols(c_idx,:));
            end
        end
        for c_idx = 1:vis_num
            plot(d.dstats.vis(c_idx).hs, d.mstats.vis(c_idx).rt_corr, ...
                '-','LineWidth',1.5,'Color',vis_cols(c_idx,:));
            plot_rt_corr(d.dstats.vis(c_idx),vis_cols(c_idx,:));
        end
        % formatting
        set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
            'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
            'XTick',[-10,0,10],'XTickLabel',{});
        ylabel('reaction time');
        
        subplot(2,1,2);  hold on;
        xlim([-17 17]);  ylim([0 1]);
        for c_idx = 1:vis_num
            plot(d.dstats.vis(c_idx).hs, d.mstats.vis(c_idx).pr, ...
                '-','LineWidth',1.5,'Color',vis_cols(c_idx,:));
            plot_pr(d.dstats.vis(c_idx),vis_cols(c_idx,:));
        end
        % formatting
        set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1],...
            'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1),...
            'XTick',[-10,0,10],'XTickLabel',{},...
            'YTick',0:025:1,'YTickLabel',{'0','','0.5','','1'});
        xlabel('heading');  ylabel('p(right)');
    end
    
    % vestibular
    if vest_num > 0
        figure;
        subplot(2,1,1);  hold on;
        if plot_incorr_rts
            plot(d.dstats.vest.hs, d.mstats.vest.rt_incorr,...
                '--','LineWidth',0.5,'Color',vest_col);
            plot_rt_incorr(d.dstats.vest,vest_col);
        end
        plot(d.dstats.vest.hs, d.mstats.vest.rt_corr, ...
            '-','LineWidth',1.5,'Color',vest_col);
        plot_rt_corr(d.dstats.vest,vest_col);
        % formatting
        set(gca,'Layer','top','TickDir','out','XTickLabel',{},...
            'FontName','Helvetica','FontSize',12);
        ylabel('reaction time');
        
        subplot(2,1,2);  hold on;
        plot(d.dstats.vest.hs, d.mstats.vest.pr, ...
            '-','LineWidth',1.5,'Color',vest_col);
        plot_pr(d.dstats.vest,vest_col);
        % formatting
        plot([0 0], [0 1], 'k--');
        plot(get(gca,'XLim'),[0.5 0.5], 'k--');
        set(gca,'Layer','top','TickDir','out',...
            'FontName','Helvetica','FontSize',12);
        ylim([0 1]);
        xlabel('heading');  ylabel('p(right)');        
    end
    
    % combined
    if comb_num > 0
        figure;
        subplot(2,1,1);  hold on;
        if plot_incorr_rts
            for c_idx = 1:comb_num
                plot(d.dstats.comb(c_idx).hs, d.mstats.comb(c_idx).rt_incorr,...
                    '--','LineWidth',0.5,'Color',comb_cols(c_idx,:));
                plot_rt_incorr(d.dstats.comb(c_idx),comb_cols(c_idx,:));
            end
        end
        for c_idx = 1:comb_num
            plot(d.dstats.comb(c_idx).hs, d.mstats.comb(c_idx).rt_corr, ...
                '-','LineWidth',1.5,'Color',comb_cols(c_idx,:));
            plot_rt_corr(d.dstats.comb(c_idx),comb_cols(c_idx,:));
        end
        % formatting
        set(gca,'Layer','top','TickDir','out','XTickLabel',{},...
            'FontName','Helvetica','FontSize',12);
        ylabel('reaction time');
        
        subplot(2,1,2);  hold on;
        for c_idx = 1:comb_num
            plot(d.dstats.comb(c_idx).hs, d.mstats.comb(c_idx).pr, ...
                '-','LineWidth',1.5,'Color',comb_cols(c_idx,:));
            plot_pr(d.dstats.comb(c_idx),comb_cols(c_idx,:));
        end
        % formatting
        plot([0 0], [0 1], 'k--');
        plot(get(gca,'XLim'),[0.5 0.5], 'k--');
        set(gca,'Layer','top','TickDir','out',...
            'FontName','Helvetica','FontSize',12);
        ylim([0 1]);
        xlabel('heading');  ylabel('p(right)');        
    end
end


%% plot everything
%% TODO



function plot_rt_incorr(dstats, plot_col)
%% function plotting the incorrect RT data for the given condition
%
% Note that the error bars are 2 SEM, rather than only one, to cover
% roughly 95% of the mass

hs_data = dstats.n_incorr >= 2;
for h_idx = 1:length(dstats.hs)
    if hs_data(h_idx)
        plot(dstats.hs(h_idx)*ones(1,2), ...
            dstats.rt_incorr(h_idx)+...
            [-1 1]*2*sqrt(dstats.rt_incorr_var(h_idx)/dstats.n_incorr(h_idx)),...
            '-','LineWidth',0.5,'Color',plot_col);
    end
end
plot(dstats.hs(hs_data),dstats.rt_incorr(hs_data),...
    's','MarkerSize',4,'LineWidth',0.5,'Color',plot_col,...
    'MarkerFaceColor',[1 1 1]);


function plot_rt_corr(dstats, plot_col)
%% function plotting the correct RT data for the given condition
%
% Note that the error bars are 2 SEM, rather than only one, to cover
% roughly 95% of the mass

hs_data = dstats.n_corr >= 2;
for h_idx = 1:length(dstats.hs)
    if hs_data(h_idx)
        plot(dstats.hs(h_idx)*ones(1,2), ...
            dstats.rt_corr(h_idx)+...
            [-1 1]*2*sqrt(dstats.rt_corr_var(h_idx)/dstats.n_corr(h_idx)),...
            '-','LineWidth',0.5,'Color',plot_col);
    end
end
plot(dstats.hs(hs_data),dstats.rt_corr(hs_data),...
    'o','MarkerSize',4,'LineWidth',0.5,'Color',plot_col,...
    'MarkerFaceColor',[1 1 1]);


function plot_pr(dstats, plot_col)
%% function plotting p(right) for given condition

if isfield(dstats,'pr_u') && isfield(dstats,'pr_l')
    for i = 1:length(dstats.hs)
        plot(dstats.hs(i)*[1 1],[dstats.pr_u(i) dstats.pr_l(i)],...
            '-','Color',plot_col);
    end
end
plot(dstats.hs,dstats.pr,'o','MarkerSize',4,'LineWidth',0.5,...
    'Color',plot_col,'MarkerFaceColor',[1 1 1]);
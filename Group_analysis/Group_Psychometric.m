function function_handles = Group_Psychometric(XlsData,if_tolerance)

%% Get data
monkey = get(findall(gcbf,'tag','monkey'),'value');

num = XlsData.num;
txt = XlsData.txt;
raw = XlsData.raw;
header = XlsData.header;

% == Figure default
set(0, 'defaultFigureRenderer', 'painters')
set(0, 'defaultFigureRendererMode', 'Manual')
set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46]/255);

colors = mat2cell([41 89 204; 248 28 83; 14 153 46]/255, ones(3,1));
heading_colors = [0 0 0; 0.3137 0.1961 0.1249; 0.6274 0.3921 0.2497; 0.9363 0.5851 0.3726; 1.0000 0.7812 0.4975];
stimtype_name = {'Vestibular'; 'Visual'; 'Combined'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======= Change the file directory===========

% ==== Get data ====
if monkey == 1 % Duncan
    
    mat_address = 'D:\Paper_rawdata\Raw_data\CN\Behavior\CN_m10&13&15_Behavior\m15\';

    mask_all = (strcmp(txt(:,header.Protocol),'BHD') | strcmp(txt(:,header.Protocol),'HD')) & (num(:,header.Monkey) == 15) & strcmp(txt(:,header.Area),'CD');
    
    
elseif monkey == 2 % Fara
    mat_address = 'D:\Paper_rawdata\Raw_data\CN\Behavior\CN_m10&13&15_Behavior\m13\';
    mask_all = (strcmp(txt(:,header.Protocol),'BHD') | strcmp(txt(:,header.Protocol),'HD')) & (num(:,header.Monkey) == 13);
    
elseif monkey == 3 % Messi
    mat_address = 'D:\Paper_rawdata\Raw_data\CN\Behavior\CN_m10&13&15_Behavior\m10_RT\';
    mask_all = (strcmp(txt(:,header.Protocol),'BHD')| strcmp(txt(:,header.Protocol),'HD') & (num(:,header.Monkey) == 10));
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xls_num = num(mask_all,:);
xls_txt = txt(mask_all,:);
xls_raw = raw(mask_all,:);

% Exclude duplications (more than one cells for one psychometric session)
[~,unique_session]=unique(xls_txt(:,header.FileNo),'stable');

xls_num = xls_num(unique_session,:);
xls_txt = xls_txt(unique_session,:);
xls_raw = xls_raw(unique_session,:);

cd(mat_address);

% group_result(sum(mask_all)).raw = [];

nSession = size(xls_raw,1);
real_Session = [];

% Load files
for i = 1:nSession
    % Find .mat file
    try
        fileName = dir([mat_address xls_txt{i,header.FileNo} '_*.mat']);
        result = load([mat_address fileName(1).name]);  % Only the first appearance
        result = result.result;
        
        if result.repetitionN >= 8 && length(result.unique_motion_coherence)==1 && (~isnan(result.Thresh_psy{3})||result.Thresh_psy{3} > 0.5) % Only repN > 8 and only have one coherence (and exclude abnormal thresholds)
            
            % Filling in group_result
            real_Session = [real_Session i];
            group_result(length(real_Session)) = struct(result);
        end
        
        
    catch
        disp(['Read data error:  ' mat_address xls_txt{i,header.FileNo}]);
        %         keyboard;
    end
    
end

nSession = length(real_Session);
xls_num = xls_num(real_Session,:);
xls_txt = xls_txt(real_Session,:);
xls_raw = xls_raw(real_Session,:);


% Some more preprocessings
for i = 1:nSession
    group_result(i).staircase = range(group_result(i).fit_data_psycho_cum{group_result(i).unique_stim_type(1)}(:,end)) > 3; % Staircase marker
end

% Pack data
staircase = [group_result.staircase]';
CR_session =  reshape(cell2mat([group_result.correct_rate]),3,[])';

if if_tolerance
    Thresh_session = reshape(cell2mat([group_result.Thresh_psy_tol]),3,[])';
    Bias_session = reshape(cell2mat([group_result.Bias_psy_tol]),3,[])';
    Thresh_shift_session = reshape([group_result.psy_thresh_shift_tol],3,[])';
else
    Thresh_session = reshape(cell2mat([group_result.Thresh_psy]),3,[])';
    Bias_session = reshape(cell2mat([group_result.Bias_psy]),3,[])';
    Thresh_shift_session = reshape([group_result.psy_thresh_shift],3,[])';
end

Coh_session = nan(nSession,1);
for i = 1:nSession
    if ~isempty([find(group_result(i).unique_stim_type == 2) find(group_result(i).unique_stim_type == 3)])
        Coh_session(i,1) = mode(group_result(i).raw(:,2));
    end
end


% Auxillary
tf_session = xls_num(:, header.HD_TargFirst);
glass_session = xls_num(:, header.Psy_glass);
platform_session = xls_num(:,header.Psy_platform);

%% Turn sessions into days (Average every day)
day = xls_num(:,1);
unique_day = unique(day);

nDay = length(unique_day);
Thresh_day = nan(size(unique_day,1),3);
Bias_day = nan(size(unique_day,1),3);
Coh_day = nan(size(unique_day,1),1);

prediction_ratio_interaction = nan(size(unique_day,1),2); % 1:interleaved; 2: comb_alone

for dd = 1:length(unique_day)
    
    thresh_today = Thresh_session(day == unique_day(dd),:);
    bias_today = Bias_session(day == unique_day(dd),:);
    coh_today = Coh_session(day == unique_day(dd),:);
    
    for condition = 1:3
        % Day-average
        Thresh_day(dd,condition) = mean(thresh_today(~isnan(thresh_today(:,condition)),condition));
        Bias_day(dd,condition) = mean(bias_today(~isnan(thresh_today(:,condition)),condition));
        
        % Auxiliary
        tf_day (dd) = mode(tf_session(day == unique_day(dd)));
        glass_day (dd) = mode(glass_session(day == unique_day(dd)));
        platform_day (dd) = mode(platform_session(day == unique_day(dd)));
    end
    
    % For interaction analysis
    sigma_pred_today = 1./sqrt(1./Thresh_day(dd,1).^2 + 1./Thresh_day(dd,2).^2);
    comb_alone_today = isnan(thresh_today(:,1)) & isnan(thresh_today(:,2)) & ~isnan(thresh_today(:,3));
    interleave_today = ~any(isnan(thresh_today),2);
    if ~isempty(interleave_today)
        prediction_ratio_interaction(dd,1) = mean(thresh_today(interleave_today,3))/sigma_pred_today;
    end
    if ~isempty(comb_alone_today)
        prediction_ratio_interaction(dd,2) = mean(thresh_today(comb_alone_today,3))/sigma_pred_today;
    end
    
    % Find coherence for each day
    Coh_day(dd,1) = mean(coh_today(~isnan(coh_today(:,1)),1));
end

interaction_pairs = fliplr(prediction_ratio_interaction(~any(isnan(prediction_ratio_interaction),2),:));

showPlatform = 0; % ~isempty(strfind(mat_address,'Polo'));
platforms = [301 109 103 102];
platformColor = {'k','r','g','g'};

%% ====================================== Function Handles =============================================%

function_handles = {
    'Training days', {
    'Threshold and bias', @f1p1;
    'Prediction ratio', @f1p2;
    'Average threshold & bias', @f1p3;
    };
    
    'Training sessions', {
    'Threshold and bias', @f2p1;
    'Prediction ratio', @f2p2;
    'Average threshold & bias', @f2p3;
    'Average reaction time', @f2p4;
    };
    
    };

%% ====================================== Function Definitions =============================================%

    function f1p1(debug)       % Threshold and Bias: Training day
        if debug;  dbstack;  keyboard;  end
        
        set(figure(76),'Pos',[10 258 880 705]); clf; hold on;
        
        % Plot threshold
        for condition = 1:3
            h(condition) = plot(1:nDay,Thresh_day(:,condition),[colors{condition} 'o'],'markerfacecolor',colors{condition});
            
            if ~all(isnan(Thresh_day(:,condition)))
                xx = find(~isnan(Thresh_day(:,condition)),1,'first'):find(~isnan(Thresh_day(:,condition)),1,'last');
                yy = smooth(Thresh_day(:,condition),30,'rloess');
                plot(xx,yy(xx),['-' colors{condition}],'linewid',2);
            end
            
            %      axis([0 50 2 190]);
        end
        
        % Plot coherence
        h(4) = plot(1:nDay,Coh_day,'rs');
        
        legend(h,{'Vestibular','Visual','Combined','Coherence'});
        set(gca,'yscale','log','ytick',[1 5 10]);
        
        xlabel('Training Day');
        ylabel('Threshold'); xlim([0 nDay+2]);
        gg = gca;
        GrayGrid(gg);
        Annotation(1:nDay,tf_day,glass_day);
        SetFigure(20);
        
        %% -- Abs(bias) --
        set(figure(77),'Pos',[67 182 887 719]); clf; hold on;
        
        set(gca,'yscale','log','ytick',[0.1 1 10],'yticklabel',[0.1 1 10]);
        %         end
        % commented by ZZ 20201205
        
        for condition = 1:3
            plot(1:nDay,abs(Bias_day(:,condition)),[colors{condition} 's'],'markerfacecolor',colors{condition});
            
            if ~all(isnan(Bias_day(:,condition)))
                xx = find(~isnan(Bias_day(:,condition)),1,'first'):find(~isnan(Bias_day(:,condition)),1,'last');
                yy = smooth(abs(Bias_day(:,condition)),30,'rloess');
                plot(xx,yy(xx),['-' colors{condition}],'linewid',2);
            end
        end
        
        xlim([0 nDay+2]);
        xlabel('Training Day');
        ylabel('|Bias|');
        
        GrayGrid(gca);
        
        SetFigure(20);
        
        %% == u time course ==
        
        set(figure(78),'position',[119 128 880 705]); clf
        plot([0 nDay],[0 0],'k--','linew',2);
        hold on
        
        for condition = 1:3
            plot(Bias_day(:,condition),[colors{condition} 'o'],'markerfacecol',colors{condition},'markersize',10);
            
            if ~all(isnan(Bias_day(:,condition)))
                xx = find(~isnan(Bias_day(:,condition)),1,'first'):find(~isnan(Bias_day(:,condition)),1,'last');
                yy = smooth(Bias_day(:,condition),30,'rloess');
                plot(xx,yy(xx),['-' colors{condition}],'linewid',2);
            end
        end
        % plot(Bias_pred_session,'go','markersize',13,'linewidt',1.5);
        xlabel('Training Day');
        ylabel('Bias');
        xlim([0 nDay+2]);
        ylim(1.1*[-max(abs(ylim)) max(abs(ylim))]);
        
        Annotation(1:nDay,tf_day,glass_day);
        
        SetFigure();
    end
    function f1p2(debug)       % Prediction ratio : Training Day
        if debug;  dbstack;  keyboard;  end
        
        sigma_pred_day = 1./sqrt(1./Thresh_day(:,1).^2 + 1./Thresh_day(:,2).^2);
        sigma_pred_ratio_day = Thresh_day(:,3)./sigma_pred_day;
        
        set(figure(79),'position',[19 85 1388 535]);
        clf;
        plot([0 nDay],[1 1],'k--','linew',2); hold on;
        
        % Min threshold
        % plot(min(Thresh_session(:,1:2),[],2)./sigma_pred,'v','color',[0.5 0.5 0.5],'markerfacecol',[0.5 0.5 0.5]);
        % plot(max(Thresh_session(:,1:2),[],2)./sigma_pred,'k^','markerfacecol','k');
        
        if ~all(isnan(sigma_pred_day))
            plot(smooth(min(Thresh_day(:,1:2),[],2)./sigma_pred_day,20,'rloess'),['-.k'],'linewid',2);
            plot(smooth(max(Thresh_day(:,1:2),[],2)./sigma_pred_day,20,'rloess'),['-.k'],'linewid',2);
        end
        
        
        % plot(1:failureBegin-1,sigma_pred_ratio(1:failureBegin-1),'ks','markerfacecol','k','markersize',9);
        % plot(failureBegin:n, sigma_pred_ratio(failureBegin:end),'ks','markersize',9);
        
        if ~showPlatform
            scatter(1:length(sigma_pred_ratio_day),sigma_pred_ratio_day,150,linspace(0,1,length(sigma_pred_ratio_day)),'fill','s');
        else
            for pp = 1:length(platforms)
                ind = find(platform_day == platforms(pp));
                plot(ind,sigma_pred_ratio_day(ind),'s','color',platformColor{pp},'markerfacecol',platformColor{pp},'markersize',9);
            end
        end
        
        if ~all(isnan(sigma_pred_ratio_day))
            plot(smooth(sigma_pred_ratio_day,30,'rloess'),['-k'],'linewid',2);
        end
        
        ylim([0.3 3]); xlim([0 nDay+2]);
        set(gca,'YScale','log','Ytick',[0.5 1 2:3]);
        
        xlabel('Training day');
        ylabel('Prediction ratio');
        % ylabel('Prediction ratio = actual / predicted threshold');
        
        Annotation(1:nDay,tf_day,glass_day);
        
        SetFigure();
    end

    function f1p3(debug)       %  Averaged threshold & bias
        % bias analysis was added by ZZ 20210203
        if debug;  dbstack;  keyboard;  end
        
        set(figure(20),'position',[19 85 1388 535]); clf;
        
        if monkey == 1
            average_duration = 57: nDay;
        elseif monkey == 2
            average_duration = 1:nDay;
        end
        
        % -- Threshold
        thres_to_average = Thresh_day(average_duration,:);
        thres_to_average(any(isnan(thres_to_average),2),:) = [];
        
        % Predicted threshold
        thres_to_average(:,4) = sqrt(thres_to_average(:,1).^2 .* thres_to_average(:,2).^2 ./ (thres_to_average(:,1).^2 + thres_to_average(:,2).^2));
        
        mean_threshold = mean(thres_to_average);
        sem_threshold = std(thres_to_average)/sqrt(size(thres_to_average,1));
        
        % -- Bias
        bias_to_average = Bias_day(average_duration, :);
        bias_to_average(any(isnan(bias_to_average),2),:) = [];
        
        mean_bias = mean(abs(bias_to_average));   % absolute value
        sem_bias = std(abs(bias_to_average)) / sqrt(size(bias_to_average, 1));
        
        
        % Plotting
        
        color_bar = {'b','r','g','c'};
        
        % -- plotting average threshold
        subplot(1,2,1);
        xlim([0.5 4.5]);
        hold on;
        for i = 1:4
            bar(i,mean_threshold(i),0.7,'facecol',color_bar{i},'edgecol','none');
            h = errorbar(i,mean_threshold(i),sem_threshold(i),color_bar{i},'linestyle','none','linewidth',3);
            errorbar_tick(h,13);
        end
        
        % Statistics
        [~,p_vest_vis] = ttest(thres_to_average(:,1)-thres_to_average(:,2));
        [~,p_comb_vest] = ttest(thres_to_average(:,1)-thres_to_average(:,3));
        [~,p_comb_vis] = ttest(thres_to_average(:,2)-thres_to_average(:,3));
        [~,p_comb_pred] = ttest(thres_to_average(:,3)-thres_to_average(:,4));
        
        set(gca,'xtick',1:4,'xticklabel',{'Vestibular','Visual','Combined','Optimal'});
        rotateXLabels(gca,45);
        ylabel('Threshold');
        text(2.5,max(ylim),sprintf('p vest-vis = %g, p comb-vest = %g\np comb-vis = %g, p comb-pred = %g',p_vest_vis,p_comb_vest,p_comb_vis,p_comb_pred),'fontsize',3);
        
        % -- plotting average bias
        % added by ZZ 20210203
        subplot(1,2,2);
        xlim([0.5 3.5]);
        hold on;
        for i = 1:3
            bar(i, mean_bias(i),0.7,'facecol',color_bar{i},'edgecol','none');
            h = errorbar(i,mean_bias(i),sem_bias(i),color_bar{i},'linestyle','none','linewidth',3);
            errorbar_tick(h,13);
        end
        
        % Statistics
        [~,p_vest_0] = ttest(bias_to_average(:,1));
        [~,p_vis_0] = ttest(bias_to_average(:,2));
        [~,p_comb_0] = ttest(bias_to_average(:,3));
        [~,p_vest_vis] = ttest(abs(bias_to_average(:,1))-abs(bias_to_average(:,2)));
        [~,p_vest_comb] = ttest(abs(bias_to_average(:,1))-abs(bias_to_average(:,3)));
        [~,p_vis_comb] = ttest(abs(bias_to_average(:,2))-abs(bias_to_average(:,3)));
        
        set(gca,'xtick',1:3,'xticklabel',{'Vestibular','Visual','Combined'});
        rotateXLabels(gca,45);
        ylabel('|Bias|');
        text(1,max(ylim),sprintf('p vest = %g, p vis = %g, p comb = %g\np vest-vis = %g, p vest-comb = %g, p vis-comb = %g',...
            p_vest_0,p_vis_0,p_comb_0,p_vest_vis,p_vest_comb,p_vis_comb),'fontsize',3);
        
        set(gcf,'defaultaxesfontsize',21);
        suptitle(sprintf('n = %g days',length(thres_to_average)));
        SetFigure(15);
    end

    function f2p1(debug)       %  Threshold and Bias: Session
        if debug;  dbstack;  keyboard;  end
        
        % == sigma time course ==
        
        set(figure(3),'position',[19 85 1388 535]); clf
        
        comb_alone =  isnan(Thresh_session(:,1)) & isnan(Thresh_session(:,2)) & ~isnan(Thresh_session(:,3)); % Combine alone vs interleaved
        
        for condition = 1:3
            plot(Thresh_session(:,condition),[colors{condition} 'o'],'markerfacecol',colors{condition},'markersize',10); hold on;
            if ~all(isnan(Thresh_session(:,condition)))
                plot(smooth(Thresh_session(:,condition),30,'rloess'),['-' colors{condition}],'linewid',2);
            end
        end
        
        % Replot comb_alone
        plot(find(comb_alone),Thresh_session(comb_alone,3),'ko','markerfacecol','g','markersize',10,'linewi',2);
        
        % plot(sigma_pred_session,'go','markersize',10,'linewidt',1.5);
        
        % plot(coh);
        xlabel('Session');
        ylabel('Threshold');
        xlim([0 nSession+2]);
        set(gca,'YScale','log','Ytick',[1:10]);
        
        Annotation(day,tf_session,glass_session);
        
        SetFigure();
        
        % == u time course ==
        
        Bias_pred_session = (Thresh_session(:,2).^2 .* Bias_session(:,1) + Thresh_session(:,1).^2 .* Bias_session(:,2))./(Thresh_session(:,1).^2 + Thresh_session(:,2).^2);
        
        set(figure(4),'position',[19 85 1388 535]); clf
        plot([0 nSession],[0 0],'k--','linew',2);
        hold on
        
        for condition = 1:3
            plot(Bias_session(:,condition),[colors{condition} 'o'],'markerfacecol',colors{condition},'markersize',10);
            if ~all(isnan(Bias_session(:,condition)))
                plot(smooth(Bias_session(:,condition),30,'rloess'),['-' colors{condition}],'linewid',2);
            end
        end
        
        % Replot comb_alone
        plot(find(comb_alone),Bias_session(comb_alone,3),'ko','markerfacecol','g','markersize',10,'linewi',2);
        
        % plot(Bias_pred_session,'go','markersize',13,'linewidt',1.5);
        xlabel('Session');
        ylabel('|Bias|');
        ylim(1.1*[-max(abs(ylim)) max(abs(ylim))]);
        xlim([0 nSession+2]);
        
        Annotation(day,tf_session,glass_session);
        
        
        SetFigure();
    end
    function f2p2(debug)       % Prediction ratio: Session
        if debug;  dbstack;  keyboard;  end
        
        sigma_pred_session = 1./sqrt(1./Thresh_session (:,1).^2 + 1./Thresh_session(:,2).^2);
        sigma_pred_ratio_session = Thresh_session(:,3)./sigma_pred_session;
        
        set(figure(5),'position',[19 85 1388 535]);
        clf;
        plot([0 nSession],[1 1],'k--','linew',2); hold on;
        
        % Min threshold
        % plot(min(Thresh_session(:,1:2),[],2)./sigma_pred,'v','color',[0.5 0.5 0.5],'markerfacecol',[0.5 0.5 0.5]);
        % plot(max(Thresh_session(:,1:2),[],2)./sigma_pred,'k^','markerfacecol','k');
        
        if ~all(isnan(sigma_pred_session))
            plot(smooth(min(Thresh_session(:,1:2),[],2)./sigma_pred_session,20,'rloess'),['-.k'],'linewid',2);
            plot(smooth(max(Thresh_session(:,1:2),[],2)./sigma_pred_session,20,'rloess'),['-.k'],'linewid',2);
        end
        
        % plot(1:failureBegin-1,sigma_pred_ratio(1:failureBegin-1),'ks','markerfacecol','k','markersize',9);
        % plot(failureBegin:n, sigma_pred_ratio(failureBegin:end),'ks','markersize',9);
        
        if ~showPlatform
            scatter(1:nSession,sigma_pred_ratio_session,150,linspace(0,1,nSession),'fill','s');
        else
            for pp = 1:length(platforms)
                ind = find(platform_session == platforms(pp));
                plot(ind,sigma_pred_ratio_session(ind),'s','color',platformColor{pp},'markerfacecol',platformColor{pp},'markersize',9);
            end
        end
        
        if ~all(isnan(sigma_pred_ratio_session))
            plot(smooth(sigma_pred_ratio_session,30,'rloess'),['-k'],'linewid',2);
        end
        
        ylim([0.3 3]); xlim([0 nSession+2]);
        set(gca,'YScale','log','Ytick',[0.5 1 2:3]);
        
        xlabel('Session');
        ylabel('Prediction ratio');
        % ylabel('Prediction ratio = actual / predicted threshold');
        
        Annotation(day,tf_session,glass_session);
        
        SetFigure();
    end

    function f2p3(debug)       %  Averaged threshold
        if debug;  dbstack;  keyboard;  end
        
        set(figure(21),'position',[19 85 1388 535]); clf;
        
        %         if monkey == 1 % Duncan
        average_duration = 1:nSession;
        %         end
        
        % -- Threshold
        thres_to_average = Thresh_session(average_duration,:);
        thres_to_average(any(isnan(thres_to_average),2),:) = [];
        
        % Predicted threshold
        thres_to_average(:,4) = sqrt(thres_to_average(:,1).^2 .* thres_to_average(:,2).^2 ./ (thres_to_average(:,1).^2 + thres_to_average(:,2).^2));
        
        mean_threshold = mean(thres_to_average);
        sem_threshold = std(thres_to_average)/sqrt(size(thres_to_average,1));
        
        % -- Bias
        bias_to_average = Bias_session(average_duration,:);
        bias_to_average(any(isnan(bias_to_average),2),:) = [];
        
        mean_bias = mean(abs(bias_to_average));   % absolute value
        sem_bias = std(abs(bias_to_average)) / sqrt(size(bias_to_average, 1));
        
        color_bar = [colors; 0 0 0]; 

        
        % -- plotting average threshold
        subplot(1,2,1);
        xlim([0.5 4.5]);
        hold on;
        for i = 1:4
            bar(i,mean_threshold(i),0.7,'facecol',color_bar{i},'edgecol','none');
            h = errorbar(i,mean_threshold(i),sem_threshold(i),'color',color_bar{i},'linestyle','none','linewidth',3);
            errorbar_tick(h,13);
        end
        
        % Statistics
        [~,p_vest_vis] = ttest(thres_to_average(:,1)-thres_to_average(:,2));
        [~,p_comb_vest] = ttest(thres_to_average(:,1)-thres_to_average(:,3));
        [~,p_comb_vis] = ttest(thres_to_average(:,2)-thres_to_average(:,3));
        [~,p_comb_pred] = ttest(thres_to_average(:,3)-thres_to_average(:,4));
        
        set(gca,'xtick',1:4,'xticklabel',{'Vestibular','Visual','Combined','Optimal'});
        rotateXLabels(gca,45);
        ylabel('Threshold');
        text(2.5,max(ylim),sprintf('p vest-vis = %g, p comb-vest = %g\np comb-vis = %g, p comb-pred = %g',p_vest_vis,p_comb_vest,p_comb_vis,p_comb_pred),'fontsize',3);
        
        % -- plotting average bias
        % added by ZZ 20210203
        subplot(1,2,2);
        xlim([0.5 3.5]);
        hold on;
        for i = 1:3
            bar(i, mean_bias(i),0.7,'facecol',color_bar{i},'edgecol','none');
            h = errorbar(i,mean_bias(i),sem_bias(i),'color',color_bar{i},'linestyle','none','linewidth',3);
            errorbar_tick(h,13);
        end
        
        % Statistics
        [~,p_vest_0] = ttest(bias_to_average(:,1));
        [~,p_vis_0] = ttest(bias_to_average(:,2));
        [~,p_comb_0] = ttest(bias_to_average(:,3));
        [~,p_vest_vis] = ttest(abs(bias_to_average(:,1))-abs(bias_to_average(:,2)));
        [~,p_vest_comb] = ttest(abs(bias_to_average(:,1))-abs(bias_to_average(:,3)));
        [~,p_vis_comb] = ttest(abs(bias_to_average(:,2))-abs(bias_to_average(:,3)));
        
        set(gca,'xtick',1:3,'xticklabel',{'Vestibular','Visual','Combined'});
        rotateXLabels(gca,45);
        ylabel('|Bias|');
        text(1,max(ylim),sprintf('p vest = %g, p vis = %g, p comb = %g\np vest-vis = %g, p vest-comb = %g, p vis-comb = %g',...
            p_vest_0,p_vis_0,p_comb_0,p_vest_vis,p_vest_comb,p_vis_comb),'fontsize',3);
        
        set(gcf,'defaultaxesfontsize',21);
        suptitle(sprintf('n = %g sessions',length(thres_to_average)));
        SetFigure(15);
    end


% Reaction time average
    function f2p4(debug)
        if debug; dbstack; keyboard; end
        
        if monkey ~= 3
            warning('Analysis of reaction time is reserved for Messi');
            return;
        end
        
        RT_headings = [-8 -4 -2 -1 0 0 1 2 4 8];   % two zero headings for left and right choice
        unique_headings =  [-8 -4 -2 -1 0 1 2 4 8];
        
        % Initiation for constructing a matrix containing all correct trial
        % For following Model Fitting
        conca_RT = [];  conca_session = [];
%         conca_outcome = []; conca_condition = []; conca_heading = []; conca_choice = [];  
        
        for h = 1:length(RT_headings)
            temp_RT_correct_median = [];
            
            for sess = 1 : size(group_result,2)
                
                if length(group_result(sess).unique_stim_type) < 3
                    continue;
                end
                
                if ismember(RT_headings(h), group_result(sess).RT_heading)   % no 0 degree
                    
                    hh = find(group_result(sess).RT_heading == RT_headings(h));    % heading location in each session
                    temp_RT_correct_median = [temp_RT_correct_median group_result(sess).RT_correct_median(:,hh)];
                    
                else
                    continue;
                end
                % Make full use of previous save data
                % Concatating all trials for DDM fitting 
                % Stim. type, Coherence, Heading, Choice, Outcome, RT
                conca_RT = [conca_RT; group_result(sess).raw];
                conca_session = [conca_session; ones(length(group_result(sess).raw),1)*sess];
                
            end
            
            if h == 5   % 0 heading, left choice
                RT_correct_median{h} = temp_RT_correct_median(:,1:2:end);
            elseif h == 6 % 0 heading, left choice
                RT_correct_median{h} = temp_RT_correct_median(:,2:2:end);
            else
                RT_correct_median{h} = temp_RT_correct_median;
            end
        end
        
        % mean of median
        mean_RT_correct_median = cellfun(@(x) mean(x,2), RT_correct_median, 'UniformOutput', 0);
        mean_RT_correct_median = cell2mat(mean_RT_correct_median);
                
        set(figure(717), 'name', 'Mean of RT median', 'pos', [50 50 1000 600]); clf;
        for condi = 1:3
            subplot(1,3,condi); hold on;
            for h = 1:length(RT_headings)
                plot(RT_headings(h)*ones(size(RT_correct_median{h},2)), RT_correct_median{h}(condi,:), 'ko','MarkerSize',6, 'MarkerFaceColor', 'none');
            end
            
            plot(RT_headings, mean_RT_correct_median(condi,:), 'o', 'color', colors{condi,:}, 'MarkerSize',9,'MarkerFaceColor',colors{condi,:});
            plot(RT_headings(1:5), mean_RT_correct_median(condi,1:5), '-', 'color', colors{condi,:}, 'LineWidth',2);
            plot(RT_headings(6:end), mean_RT_correct_median(condi,6:end), '-', 'color', colors{condi,:}, 'LineWidth',2);
            
            xlim([min(RT_headings)-1 max(RT_headings)+1]);
            xticks(unique_headings);
            ylim([650 1500]);
            
            if condi == 1
                ylabel('Reaction time (ms)');
            elseif condi == 2
                title([num2str(length(RT_correct_median{h})) ' Sessions'])
                xlabel('Heading (\circ)')
            end
        end
        SetFigure();
        
        
        % ===============  Concatenated All-Trial Information ==========================
        left_is_correct = ((conca_RT(:,4)==1) & (conca_RT(:,5)==0)) | ((conca_RT(:,4)==2) & (conca_RT(:,5)==5)); 
        
        ALL_trial_info = [conca_session, conca_RT, left_is_correct]; 
        
        cd(mat_address)
        if exist('ALL_trial_info.mat')
            delete('ALL_trial_info.mat');
        end
        % Session#; Condition#; Heading; RT; Choices; Outcomes;
        disp('..............Save all-trial RT information..............')
        save('ALL_trial_info.mat', 'ALL_trial_info');
        
        if exist('ALL_trial_info.csv')
            delete('ALL_trial_info.csv');
        end
        disp('...............Convert into csv. file .........................')
        titles = {'number','session', 'condition', 'heading', 'rt', 'choice', 'outcome', 'left_is_correct'};
        result_table = table((1:length(ALL_trial_info))', ALL_trial_info(:,1), ALL_trial_info(:,2), ALL_trial_info(:,4),...
            ALL_trial_info(:,7)/1000,ALL_trial_info(:,5), ALL_trial_info(:,6)==0, ALL_trial_info(:,8), 'VariableNames', titles);   % unit of rt is converted into s
        writetable(result_table, 'ALL_trial_info.csv');
       %========================================================================================
       
    end

%% Improvement v.s. correct rate

    function p = p_value_binomial(totalN, n, prop)
        if nargin < 3
            prop = 0.5;
        end
        
        x = 0:1:totalN;
        y = binopdf(x,totalN,prop);
        
        if n < 1  % Prop
            n = totalN * n;
        end
        
        % Putting this ahead makes the p_value symmetric (espeically when totalN is small).  HH20150125
        if n > totalN/2
            n = totalN-n;
        end
        
        p = sum(y(1:min(find(x>=n))));
        
    end

    function output = p_value_constantHD(headingN, repN, prop)
        simN = 1500;
        xcenterN = 200;
        
        unique_heading = [-fix(headingN/2):-1 , 1:fix(headingN/2)];
        
        headings = zeros(length(unique_heading),repN*simN);
        
        for i = 1:repN*simN
            headings(:,i) = unique_heading(randperm(length(unique_heading)));
        end
        
        headings = reshape(headings,length(unique_heading)*repN,[]);
        shift_prop = sum(sign(headings(1:end-1,:)) ~= sign(headings(2:end,:)))/(size(headings,1)-1);
        
        [n,x] = hist(shift_prop,xcenterN);
        
        output.mean = mean(shift_prop);
        
        for pp = 1:length(prop)
            output.p(pp) = sum(n(1:min(find(x>prop(pp)))))/ simN;
            if prop >= output.mean
                output.p = 1-output.p;
            end
        end
        
    end

    function Annotation(date,tf,glass)
        
        % Same day
        diff0 = find(diff(date)== 0);
        
        if ~isempty(diff0)
            diff0_begs = [diff0(1); diff0(1+find(diff(diff0)>1))];
            diff0_ends = [diff0(diff(diff0)>1) ;diff0(end)]+1;
        else
            diff0_begs = [];
        end
        
        ylims = ylim;
        xlims = xlim;
        
        for j = 1:length(diff0_begs);
            plot([diff0_begs(j)-0.2 diff0_ends(j)+0.2],[ylims(1) ylims(1)],'k-','linewid',10);
        end
        
        % Period that Polo refused to work
        % plot([70 70],ylims,'--k','linewi',2);
        % plot([74 74],ylims,'--k','linewi',2);
        
        
        % -- Target first v.s. 3D glass --
        a1 = gca;
        pos = get(gca,'Position');
        a2 = axes('position',[pos(1) pos(2)+pos(4) pos(3) 0.03]);
        hold on;
        
        if ~all(isnan(tf))
            plot(find(~isnan(tf) & tf ~= 0),0,'cs','markerfacecolor','c');
        end
        
        if ~all(isnan(glass)) && ~isempty(find(glass>0))
            plot(find(glass>0),1,'ms','markerfacecolor','m');
        end
        
        axis off;
        xlim(xlims); ylim([0 1])
        
        hlink = linkprop([a1 a2],'xlim');
        set(gcf,'UserData',hlink); % Store the hlink variable in an object's UserData property to keep the links! HH20150609
        
    end

end
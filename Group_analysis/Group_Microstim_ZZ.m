% Population analysis for microstimulation results
% Running this code, the current folder need contain all microstimulation data
% Just a quick look

function Group_Microstim_ZZ(~)

pathname = uigetdir(cd, 'Choose a folder');
if pathname == 0
    msgbox('You did not choose a correct folder');
    return;
else
    cd(pathname);
end

%% Pre-definition
% Condition colors
stim_type = {'Vestibular', 'Visual', 'Combined'};
% con_color = {'b', 'r', 'g'};
con_color = [41 89 204; 248 28 83; 14 153 46]/255; 
% == Figure default
set(0, 'defaultFigureRenderer', 'painters')
set(0, 'defaultFigureRendererMode', 'Manual')
set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46]/255);

unique_heading = [-8 -4 -2 -1 0 1 2 4 8];

%% Pool all sessions
fileName_list = dir('*_ustim.mat');

% % Initiation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delta_bias= nan(size(fileName_list,1), 3);  % hist and bar function ignore NAN, so this allow blocks containing less than 3 conditions
% p_bias =delta_bias;
% delta_thresh = delta_bias;
% p_thresh = delta_bias;
% bias0 = delta_bias;  bias1 = delta_bias; 
% thresh0 = delta_bias; thresh1 = delta_bias;
% 
% % Correct rate.   Added by @ZZ 20220425 
% correct_rate0 = nan(size(fileName_list,1),3,length(unique_heading));   correct_rate1 = correct_rate0;    
% delta_correctrate = correct_rate0;
% % Add across-repetition change of bias, threshold, correct rate
% % ZZ @ 20220425
% temporal_bias0 = cell(size(fileName_list,1), 3);  temporal_bias1 = temporal_bias0;
% temporal_thresh0 = temporal_bias0;  temporal_thresh1= temporal_bias0;
% temporal_correctrate0 = temporal_bias0; temporal_correctrate1 = temporal_bias0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

counter = 0;  % total file number
for sess = 1:length(fileName_list)
    
    load(fileName_list(sess).name);
    if result.repetitionN < 8
        continue;
    end
    counter = counter +1;
    
    for condi = 1:length(result.unique_stim_type)
        bias0(counter,condi) = result.Bias_psy{1,condi};   % no-microstim control
        bias1(counter,condi) = result.Bias_psy{2, condi};  % microstim condition
        thresh0(counter,condi) = result.Thresh_psy{1,condi}; 
        thresh1(counter,condi) = result.Thresh_psy{2,condi};
        correct_rate0(counter,condi,:) = result.correct_rate{1}(condi,:);
        correct_rate1(counter,condi,:) = result.correct_rate{2}(condi,:);
        delta_correctrate(counter,condi,:) = correct_rate1(counter,condi,:) - correct_rate0(counter,condi,:);   %Positive indicate better performance
        
        delta_bias(counter,condi) = bias0(counter,condi) - bias1(counter,condi);   % positive value indicates contralateral (i.e. rightward)
        p_bias(counter,condi) = result.P_bias(condi);
        
        delta_thresh(counter,condi) = thresh1(counter,condi) - thresh0(counter,condi);   % Positive indicate worse sensitivity
        p_thresh(counter,condi) = result.P_slope(condi);
        
        % Temporal changing parameters
        temporal_bias0{counter,condi} = result.psy_bias_shift{condi,1};
        temporal_bias1{counter,condi} = result.psy_bias_shift{condi,2};
        temporal_thresh0{counter,condi} = result.psy_thresh_shift{condi,1};
        temporal_thresh1{counter,condi} = result.psy_thresh_shift{condi,2};
        temporal_correctrate0{counter,condi} = result.correct_rate_shift{condi,1};
        temporal_correctrate1{counter,condi} = result.correct_rate_shift{condi,2};
    end
    
end

monkey_ind = str2num(result.FILE(2:3)); 
if monkey_ind == 13
    markers = '^';
elseif monkey_ind == 15
    markers = 'o';
end
%% Session-by-session Bias and threshold comparison

BIAS = 1; THRESHOLD = 2;
% Pool bias and threshold
psycho_params0 = [bias0 thresh0];
psycho_params1 = [bias1 thresh1]; 
p_bias; p_thresh; 

% % mean value 
% psycho_params0_mean = nanmean(psycho_params0); 
% psycho_params1_mean = nanmean(psycho_params1);

% mean value 
psycho_params0_median = nanmedian(psycho_params0); 
psycho_params1_median = nanmedian(psycho_params1);


figure(817); clf; 
set(figure(817), 'name', 'Session-by-session comparison', 'pos',[20 10 1000 1000]);
hb = tight_subplot(3,2,[.1 .05],[.1 .1], .1); 

for option = BIAS : THRESHOLD
    for condi = 1:3
        % paried t test
%         [~, pp(condi+(option-1)*3)] = ttest(psycho_params0(:,condi+(option-1)*3), psycho_params1(:,condi+(option-1)*3));
        pp(condi+(option-1)*3) = signrank(psycho_params0(:,condi+(option-1)*3), psycho_params1(:,condi+(option-1)*3));
        axes(hb(condi+(option-1)*3)); hold on; 
        
        % Mark dots with significant difference
        if option == BIAS
            sig_ind = p_bias(:,condi) < 0.05;
        else
            sig_ind = p_thresh(:,condi) < 0.05;
        end
        
        % plot
        plot(psycho_params0(:, condi+(option-1)*3), psycho_params1(:,condi+(option-1)*3),...
            markers, 'Color', con_color(condi,:), 'MarkerSize',8,'LineWidth',0.75);
        plot(psycho_params0(sig_ind, condi+(option-1)*3), psycho_params1(sig_ind,condi+(option-1)*3),...
            markers, 'Color', con_color(condi,:), 'MarkerSize',8, 'MarkerFaceColor', con_color(condi,:),'LineWidth',0.75); 
        
        % Line of unity
        min_scale = min([psycho_params0(:, condi+(option-1)*3);psycho_params1(:,condi+(option-1)*3)]);
        max_scale = max([psycho_params0(:, condi+(option-1)*3);psycho_params1(:,condi+(option-1)*3)]);
        axis([min_scale-0.2 max_scale*1.05 min_scale-0.2 max_scale*1.05]);
        plot([min_scale-0.2 max_scale*1.05], [min_scale-0.2 max_scale*1.05], 'k--');
        
        % Anotate mean value
        plot(psycho_params0_median(:, condi+(option-1)*3),min_scale+0.1, 'v', 'Color', con_color(condi,:), 'MarkerSize', 12);
        text(psycho_params0_median(:, condi+(option-1)*3),min_scale+0.1,...
            num2str(psycho_params0_median(:, condi+(option-1)*3)),...
            'Color', con_color(condi,:), 'FontSize', 15, 'FontWeight', 'bold')     % Text location is always changed after axis square !!
        
        plot(min_scale+0.1, psycho_params1_median(:, condi+(option-1)*3), '<', 'Color', con_color(condi,:), 'MarkerSize', 12);
        text(min_scale+0.1, psycho_params1_median(:, condi+(option-1)*3),...
            num2str(psycho_params1_median(:, condi+(option-1)*3)),...
            'Color', con_color(condi,:), 'FontSize', 15, 'FontWeight', 'bold')     % Text location is always changed after axis square !!
%         text(psycho_params0_mean(:, condi+(option-1)*3),min_scale+0.1, '\downarrow', 'Color', con_color(condi,:), 'FontSize', 15, 'FontWeight', 'bold');
%         text(min_scale+0.1, psycho_params1_mean(:, condi+(option-1)*3),'\leftarrow', 'Color', con_color(condi,:), 'FontSize', 15, 'FontWeight', 'bold');

        nticks = 5;
%         xticks(roundn(linspace(min_scale,max_scale,nticks),-1));
%         yticks(roundn(linspace(min_scale,max_scale,nticks),-1));
        
        txt = {[ 'n = ', num2str(sum(~isnan(p_bias(:,condi)))) ], ['p: ', num2str(roundn(pp(condi+(option-1)*3),-4))]};
        str = text(min_scale-0.2+0.05,max_scale, txt, 'HorizontalAlignment','left');
        set(str, 'Color', [pp(condi+(option-1)*3)<0.05,(pp(condi+(option-1)*3)<0.05)*0.7,0]); 
                
        % Annotation
        if condi == 1
            if option ==BIAS
                title('Bias','Color', 'k', 'FontSize',20,'FontWeight','bold');
            else
                title('Threshold','Color','k', 'FontSize',20,'FontWeight','bold');
            end
        end
        
        if condi ==2
            ylabel('\muStim (o)');
        end
        
        if condi == 3
            xlabel('No \muStim (o)');
        end
        
        axis square; 
    end
end
SetFigure();


%% Bias and threshold shift bar

% Pool delta bias and threshold
delta = [delta_bias delta_thresh];  p_value = [p_bias p_thresh];

% t-test
mean_delta = nanmean(delta);
[~, p_population] = ttest(delta);

% significant indication
Sig = p_value < 0.05;  NSig = p_value >= 0.05;

% Plotting
figure(815); clf; 
set(figure(815), 'Name',['CD_Microstim_Population' ], 'Position',[40 10 1000 1000]);
hd = tight_subplot(3,2,[.1 .05],[.1 .1], .1);

nbin = 10;
for option = BIAS : THRESHOLD
    for condi = 1: 3
        axes(hd(condi+(option-1)*3));  hold on; 
        
        [~, xcenters]  = hist(delta(:, condi+(option-1)*3), nbin);
        
        % Centralizing
        max_abs_xcenters = max(abs(xcenters));
        xcenters = linspace(-max_abs_xcenters, max_abs_xcenters, nbin);
        
        histS= hist(delta(Sig(:,condi+(option-1)*3), condi+(option-1)*3), xcenters);
        histNS = hist(delta(NSig(:, condi+(option-1)*3),condi+(option-1)*3), xcenters);
        
        hbars = bar(gca, xcenters, [histS' histNS'],1,'stacked','LineWidth',2);
        set(hbars,'EdgeColor',con_color(condi,:),'FaceColor',con_color(condi,:));
        set(hbars(2),'FaceColor','none');
                
        % Annotation
        if condi==3 && option==1
                xlabel('PSE Shift (o)', 'FontWeight', 'bold');
        end
        
        if condi==3 && option==2
            xlabel('Threshold Change (o)','FontWeight', 'bold');
        end
        
        if condi==2 && option==1
                ylabel('Number of case','FontWeight', 'bold');
        end
        
        xlim([min(xcenters), max(xcenters)]*1.2);
        ylim([0, max(histS+histNS)+1]);
        xticks(linspace(roundn(min(xcenters),-1),roundn(max(xcenters),-1),nticks));
        
        if option == 2
            txt =['n = ', num2str(sum(~isnan(delta_bias(:,condi))))];
            text(max(xcenters), max(histS+histNS), txt, 'Color', con_color(condi,:), 'HorizontalAlignment', 'right');
        end
        
        plot(mean_delta(condi+(option-1)*3)*[1 1], [0 max(histS+histNS)+1],...
            'Color', [p_population(condi+(option-1)*3)<0.05,0.7*(p_population(condi+(option-1)*3)<0.05),0],'LineWidth', 2);
        str = text(xcenters(3), max(histS+histNS)+1.5,...
            sprintf('Mean: %2.2f   p: %0.3f', mean_delta(condi+(option-1)*3), p_population(condi+(option-1)*3)));
        set(str, 'Color', [p_population(condi+(option-1)*3)<0.05,0.7*(p_population(condi+(option-1)*3)<0.05),0]); 
                
    end
end

SetFigure(15);

%% Bias and threshold comparison 
% Bar plotting @20211018
set(figure(1018),'Name', 'Bias and threshold comparison','pos',[60 10 1000 1000]);clf;
hs = tight_subplot(3,2,[.1 .05],[.05 .1],[.1 .05]);

for option = BIAS : THRESHOLD
    for cc = 1:size(bias0,2)
        min_all = min(min([psycho_params0(:,cc+3*(option-1)) psycho_params1(:,cc+3*(option-1))])); 
        max_all = max(max([psycho_params0(:,cc+3*(option-1)) psycho_params1(:,cc+3*(option-1))]));
        x_centers = linspace(min_all,max_all,nbin);
        
        hh{cc+3*(option-1)} = HistComparison({psycho_params0(:,cc+3*(option-1)) psycho_params1(:,cc+3*(option-1))},...
            'XCenters',x_centers, 'FigN',1018, 'Axes',hs(cc+3*(option-1)), 'Style','grouped');
        
        % Annotation
        if option == 1
            ylabel(stim_type{cc});
        end
        
        if cc == 3
            if option == 1
                xlabel('Bias (o)');
            else
                xlabel('Threshold (o)');
            end
        end
        
        if option == 2 && cc ==1
           lgd = legend(hs(cc+3*(option-1)), 'No \muStim', '\muStim');
        end
            
        
    end
end
SetFigure();


%% Bias and threshold change in different conditions
% t-test
figure(818); clf;
set(figure(818), 'Name',['Across-condition comparison' ], 'Position',[80 10 1000 1000]); 
ha = tight_subplot(3,2,[.1 .05],[.1 .1], .1);

comp_pair = {[1 2]; [1 3]; [2 3]};

for option = BIAS:THRESHOLD
    for cp = 1:length(comp_pair)
        axes(ha(cp+(option-1)*3)); hold on;
                
        plot(delta(:, comp_pair{cp}(1)+(option-1)*3), delta(:, comp_pair{cp}(2)+(option-1)*3),...
            ['k' markers], 'MarkerSize', 8, 'LineWidth', 0.75);
        
        
        % significance indication
        Sig1= Sig(:,comp_pair{cp}(1)+(option-1)*3) & NSig(:,comp_pair{cp}(2)+(option-1)*3);  % Only the former condition with significant difference
        Sig2 = NSig(:,comp_pair{cp}(1)+(option-1)*3) & Sig(:,comp_pair{cp}(2)+(option-1)*3);   % Only the latter
        Sig_both = Sig(:,comp_pair{cp}(1)+(option-1)*3) & Sig(:,comp_pair{cp}(2)+(option-1)*3);  % both
        
        plot(delta(Sig1,comp_pair{cp}(1)+(option-1)*3), delta(Sig1, comp_pair{cp}(2)+(option-1)*3),...
            ['k' markers],'MarkerFaceColor',con_color(comp_pair{cp}(1),:), 'MarkerSize',8, 'LineWidth',0.75);   
        plot(delta(Sig2,comp_pair{cp}(1)+(option-1)*3), delta(Sig2, comp_pair{cp}(2)+(option-1)*3),...
            ['k' markers],'MarkerFaceColor',con_color(comp_pair{cp}(2),:), 'MarkerSize',8, 'LineWidth',0.75);   
        plot(delta(Sig_both,comp_pair{cp}(1)+(option-1)*3), delta(Sig_both, comp_pair{cp}(2)+(option-1)*3),...
            ['k' markers],'MarkerFaceColor',[1,0.7,0], 'MarkerSize',8, 'LineWidth',0.75);   

        
        % Line of unity
        max_scale = max(max(delta(:, comp_pair{cp}+(option-1)*3)));
        min_scale = min(min(delta(:, comp_pair{cp}+(option-1)*3)));
        axis([min_scale-0.2 max_scale*1.05 min_scale-0.2 max_scale*1.05]);
        plot([min_scale-0.2 max_scale*1.05], [min_scale-0.2 max_scale*1.05], 'k--');
        
        plot([0 0],[min_scale-0.2 max_scale*1.05], 'k--');
        plot([min_scale-0.2 max_scale*1.05], [0 0], 'k--');
        
        nticks = 5;
        xticks(roundn(linspace(min_scale,max_scale,nticks),-1));
        yticks(roundn(linspace(min_scale,max_scale,nticks),-1));
        
        % t-test
%         [~,pp] = ttest(delta(:, comp_pair{cp}(1)+(option-1)*3), delta(:, comp_pair{cp}(2)+(option-1)*3));
        pp = signrank(delta(:, comp_pair{cp}(1)+(option-1)*3), delta(:, comp_pair{cp}(2)+(option-1)*3));
        % Linear correlation
        nan_delta = delta(all(~isnan(delta),2),:);
        [r, p] = corrcoef(nan_delta(:, comp_pair{cp}(1)+(option-1)*3), nan_delta(:, comp_pair{cp}(2)+(option-1)*3));
        [para, S] = polyfit(nan_delta(:, comp_pair{cp}(1)+(option-1)*3), nan_delta(:, comp_pair{cp}(2)+(option-1)*3), 1);
        xx = min_scale : 0.1: max_scale;
        Y = polyval(para,xx);
        hl = plot(xx,Y,'k-','LineWidth',2);


        txt_ttest = {[ 'n = ', num2str(min(sum(~isnan(delta(:,comp_pair{cp}+(option-1)*3))))) ], ['p_{ttest}: ', num2str(roundn(pp,-3))]};
        str_ttest = text(min_scale-0.2+0.05,max_scale, txt_ttest, 'HorizontalAlignment','left');
        set(str_ttest, 'Color', [pp<0.05,(pp<0.05)*0.7,0]); 
        tight_subplot(3,1,[.1 .05],[.1 .1], .1)
        txt_corr = {['r = ', num2str(roundn(r(2),-3))], ['p: ', num2str(roundn(p(2),-3))]};
        str_corr = text(max_scale-0.8, min_scale+0.5, txt_corr,'HorizontalAlignment','left');
        

        % Annotation
        if cp == 1
            if option ==BIAS
                title('\DeltaBias','Color', 'k', 'FontSize',20,'FontWeight','bold');
                ylabel('Visual','Color',con_color(comp_pair{cp}(2),:)); xlabel('Vestibular','Color',con_color(comp_pair{cp}(1),:));
            else
                title('\DeltaThreshold','Color','k', 'FontSize',20,'FontWeight','bold');
                xlabel('Vestibular','Color',con_color(comp_pair{cp}(1),:));
            end
        end
        
        if cp ==2
            if option ==BIAS
                ylabel('Combined','Color',con_color(comp_pair{cp}(2),:));
                xlabel('Vestibular','Color',con_color(comp_pair{cp}(1),:));
            else
                xlabel('Vestibular','Color',con_color(comp_pair{cp}(1),:));
            end
        end
        
        if cp == 3
            if option ==BIAS
                ylabel('Combined','Color',con_color(comp_pair{cp}(2),:));
                xlabel('Visual','Color',con_color(comp_pair{cp}(1),:));
            else
                xlabel('Visual','Color',con_color(comp_pair{cp}(1),:));
            end
        end
        
        axis square; 
        
    end
end
SetFigure();


%% Normalized Bias shift comparison in different conditions
% Added @ 20211116

% Two normalized ways
% 1. dPSE is divided by threshold in no-ustim condition of conrresponding modality within the same session (Referring to Yu & Gu, 2019, Neuron)
% 2. dPSE in each session is divided by standard deviation of all sessions
Thresh = 1;
Variance = 2; 

figure(1116); clf;
set(figure(1116), 'Name',['Across-condition Normalized dPSE Comparison' ], 'Position',[100 10 1000 1000]); 
ha = tight_subplot(3,2,[.1 .05],[.1 .1], .1);

comp_pair = {[1 2]; [1 3]; [2 3]};

for option = Thresh : Variance
    
    if option == Thresh
        norma_dPSE = delta_bias ./ thresh0;
    elseif option == Variance
        sd_dPSE = nanstd(delta_bias);
        norma_dPSE = delta_bias ./ sd_dPSE;
    end
        
    for cp = 1:length(comp_pair)
        axes(ha(cp+(option-1)*3)); hold on;
                
        plot(norma_dPSE(:, comp_pair{cp}(1)), norma_dPSE(:, comp_pair{cp}(2)),...
            ['k' markers], 'MarkerSize', 8, 'LineWidth', 0.75);
        
        % significance indication
        Sig1= Sig(:,comp_pair{cp}(1)) & NSig(:,comp_pair{cp}(2));  % Only the former condition with significant difference
        Sig2 = NSig(:,comp_pair{cp}(1)) & Sig(:,comp_pair{cp}(2));   % Only the latter
        Sig_both = Sig(:,comp_pair{cp}(1)) & Sig(:,comp_pair{cp}(2));  % both
        
        plot(norma_dPSE(Sig1,comp_pair{cp}(1)), norma_dPSE(Sig1, comp_pair{cp}(2)),...
            ['k' markers],'MarkerFaceColor',con_color(comp_pair{cp}(1),:), 'MarkerSize',8, 'LineWidth',0.75);   
        plot(norma_dPSE(Sig2,comp_pair{cp}(1)), norma_dPSE(Sig2, comp_pair{cp}(2)),...
            ['k' markers],'MarkerFaceColor',con_color(comp_pair{cp}(2),:), 'MarkerSize',8, 'LineWidth',0.75);   
        plot(norma_dPSE(Sig_both,comp_pair{cp}(1)), norma_dPSE(Sig_both, comp_pair{cp}(2)),...
            ['k' markers],'MarkerFaceColor',[1,0.7,0], 'MarkerSize',8, 'LineWidth',0.75);   

        % Line of unity
        max_scale = max(max(norma_dPSE(:, comp_pair{cp})));
        min_scale = min(min(norma_dPSE(:, comp_pair{cp})));
        axis([min_scale-0.2 max_scale*1.05 min_scale-0.2 max_scale*1.05]);
        plot([min_scale-0.2 max_scale*1.05], [min_scale-0.2 max_scale*1.05], 'k--');
        
        plot([0 0],[min_scale-0.2 max_scale*1.05], 'k--');
        plot([min_scale-0.2 max_scale*1.05], [0 0], 'k--');

        nticks = 5;
        xticks(roundn(linspace(min_scale,max_scale,nticks),-1));
        yticks(roundn(linspace(min_scale,max_scale,nticks),-1));
        
        % t-test
        [~,pp] = ttest(norma_dPSE(:, comp_pair{cp}(1)), norma_dPSE(:, comp_pair{cp}(2)));
        
        % Linear correlation
        nan_delta = norma_dPSE(all(~isnan(norma_dPSE),2),:);
        [r, p] = corrcoef(nan_delta(:, comp_pair{cp}(1)), nan_delta(:, comp_pair{cp}(2)));
        [para, S] = polyfit(nan_delta(:, comp_pair{cp}(1)), nan_delta(:, comp_pair{cp}(2)), 1);
        xx = min_scale : 0.1: max_scale;
        Y = polyval(para,xx);
        hl = plot(xx,Y,'k-','LineWidth',2);


        txt_ttest = {[ 'n = ', num2str(min(sum(~isnan(norma_dPSE(:,comp_pair{cp}))))) ], ['p_{ttest}: ', num2str(roundn(pp,-3))]};
        str_ttest = text(min_scale-0.2+0.05,max_scale, txt_ttest, 'HorizontalAlignment','left');
        set(str_ttest, 'Color', [pp<0.05,(pp<0.05)*0.7,0]); 
        
        txt_corr = {['r = ', num2str(roundn(r(2),-3))], ['p: ', num2str(roundn(p(2),-3))]};
        str_corr = text(max_scale-0.2, min_scale+0.2, txt_corr,'HorizontalAlignment','left');
        
                
        % Annotation
        if cp == 1
            if option ==Thresh
                title('Dividing Thresh','Color', 'k', 'FontSize',20,'FontWeight','bold');
                ylabel('Visual','Color',con_color(comp_pair{cp}(2),:)); xlabel('Vestibular','Color',con_color(comp_pair{cp}(1),:));
            else
                title('Dividing sd','Color','k', 'FontSize',20,'FontWeight','bold');
                xlabel('Vestibular','Color',con_color(comp_pair{cp}(1),:));
            end
        end
        
        if cp ==2
            if option ==Thresh
                ylabel('Combined','Color',con_color(comp_pair{cp}(2),:));
                xlabel('Vestibular','Color',con_color(comp_pair{cp}(1),:));
            else
                xlabel('Vestibular','Color',con_color(comp_pair{cp}(1),:));
            end
        end
        
        if cp == 3
            if option ==Thresh
                ylabel('Combined','Color',con_color(comp_pair{cp}(2),:));
                xlabel('Visual','Color',con_color(comp_pair{cp}(1),:));
            else
                xlabel('Visual','Color',con_color(comp_pair{cp}(1),:));
            end
        end
        
        axis square; 
        
    end
end
SetFigure();



%% Bias (with or without normalization) comparison across three conditions
% Added by ZZ @ 20220301, for better visualization
% Bar comparison, one-way Anova
p_dPSE = anova1(delta_bias, stim_type, 'off');
norma_dPSE = delta_bias ./ thresh0;   % Normalized by Threshold in no-microstim codition
p_norm_dPSE = anova1(norma_dPSE,stim_type,'off');

set(figure(31), 'pos', [60 10 800 900], 'Name','dPSE Comparison Across Conditions'); clf; 
hs = subplot(2,1,1); 
plot(delta_bias','-o','Color',[0.4 0.4 0.4],'LineWidth',1, 'MarkerSize',8); hold on;
for c = 1:size(delta_bias,2)
plot(c,delta_bias(:,c), 'o','LineWidth',1, 'MarkerEdgeColor',con_color(c,:), 'MarkerSize',8);
end
ylims = [min(min(delta_bias)) max(max(delta_bias))];

hs.XLim = [1-0.2 size(delta_bias,2)+0.2];
hs.YLim = [floor(ylims(1)) ceil(ylims(2))];
hs.XTick = (1:size(delta_bias,2));
hs.YTick = floor(ylims(1)):ceil(ylims(2));
hs.XTickLabel = stim_type;
ylabel('dPSE');
text(size(delta_bias,2)+0.2,ceil(ylims(2)),num2str(roundn(p_dPSE,-4)));
% Normalized dPSE by dividing threshold 
hs = subplot(2,1,2); 
plot(norma_dPSE','-o','Color',[0.4 0.4 0.4],'LineWidth',1, 'MarkerSize',8); hold on;
for c = 1:size(norma_dPSE,2)
plot(c,norma_dPSE(:,c), 'o','LineWidth',1, 'MarkerEdgeColor',con_color(c,:), 'MarkerSize',8);
end
ylims = [min(min(norma_dPSE)) max(max(norma_dPSE))];

hs.XLim = [1-0.2 size(norma_dPSE,2)+0.2];
hs.YLim = [floor(ylims(1)) ceil(ylims(2))];
hs.XTick = (1:size(norma_dPSE,2));
hs.YTick = floor(ylims(1)):ceil(ylims(2));
hs.XTickLabel = stim_type;
ylabel('Normalized dPSE');
text(size(delta_bias,2)+0.2,ceil(ylims(2)),num2str(roundn(p_norm_dPSE,-4)));

SetFigure();


%% Correct rate comparison
correctrate0 = nanmean(correct_rate0, 3);   % Avareage across different headings
correctrate1 = nanmean(correct_rate1, 3);

h = LinearCorrelation({correctrate0(:,1), correctrate0(:,2), correctrate0(:,3)},...
    {correctrate1(:,1), correctrate1(:,2), correctrate1(:,3)},...
    'Xlabel', 'Correct rate (No)' , 'Ylabel','Correct rate (\muStim)',...
    'LineStyles', {'b:','r:','g:'}, 'Markers',{markers}, 'MarkerSize',11,...
    'figN', 425,'SameScale',1, 'AxisSquare',1, 'Diagonal',1);hold on;

delete([h.group(1:3).line]);

% Anotate significant dPSE cells and mean value
correctrate0_mean = nanmean(correctrate0);   % mean correct rate of each condition
correctrate1_mean = nanmean(correctrate1); 

xlims = xlim; ylims = ylim; 
for c = 1:3
    plot(correctrate0(Sig(:,c),c), correctrate1(Sig(:,c),c), markers,'Color',con_color(c,:), 'MarkerFaceColor',con_color(c,:), 'MarkerSize',11);    
    
    plot(correctrate0_mean(c),ylims(1)+0.01, 'v','Color', con_color(c,:), 'MarkerSize', 13);
    plot(xlims(1)+0.01, correctrate1_mean(c), '<','Color', con_color(c,:), 'MarkerSize', 13);
%     text(correctrate0_mean(c),ylims(1)+0.01, '\downarrow', 'Color', con_color(c,:), 'FontSize', 20, 'FontWeight', 'bold');
%     text(xlims(1)+0.01, correctrate1_mean(c),'\leftarrow', 'Color', con_color(c,:), 'FontSize', 20, 'FontWeight', 'bold');
end

% paired t-test
[~, pp] = ttest(correctrate0, correctrate1);
text(1,1,['p: ' num2str(pp)]);

SetFigure();


% Delta correct rate in each heading
for c = 1:3
    for hh = 1:length(unique_heading)
        delta_correctrate_mean(c,hh) = nanmean(delta_correctrate(:,c,hh));
        delta_correctrate_sem(c,hh) = nanstd(delta_correctrate(:,c,hh),0,1) / sqrt(sum(~isnan(delta_correctrate(:,c,hh))));
        
        [~, p_delta_corrate(c,hh)] = ttest(delta_correctrate(:,c,hh));
    end
end

set(figure(426), 'pos', [50 60 1000, 500], 'Name', 'Delta Correct Rate'); clf; hold on;
errorbar(repmat(unique_heading',1,3), delta_correctrate_mean',delta_correctrate_sem',...
    'LineWidth', 2)
xlim([-9 9]); ylims = ylim; 
plot([-9 9], [0 0],'k--', 'linew',1);
xlabel('Heading (o)'); ylabel('\Delta Correct Rate');
legend(stim_type, 'AutoUpdate','off');

% significance indication
sig_ind = p_delta_corrate <= 0.05;  % significance index
for c = 1:3
    
    if any(sig_ind(c,:))
        text(unique_heading(sig_ind(c,:)), repmat(ylims(2)+(c-1)*0.005, 1,sum(sig_ind(c,:))),...
            '\ast', 'color', con_color(c,:), 'FontWeight','bold');
    end
end
SetFigure();
%% Temporal dynamics of deltaPSE with trials going
% Added by ZZ@ 20220425
% Time window = 3 repetitions

% Calcute mean and sem of each time windows in different conditions
for tt = 1: min(min(cellfun(@(x) size(x,2), temporal_bias0)))    
    
        bias0_temp = cellfun(@(x) x(:,tt), temporal_bias0);
        bias1_temp = cellfun(@(x) x(:,tt), temporal_bias1);
        % Delta bias
        delta_bias_temporal = bias1_temp-bias0_temp;
        delta_bias_temporal_mean(tt,:) = nanmean(delta_bias_temporal);
        delta_bias_temporal_sem(tt,:) = nanstd(delta_bias_temporal) ./ sqrt(sum(~isnan(delta_bias_temporal)));
        % For ANOVA
        delta_bias4anova(:,:,tt) = delta_bias_temporal;

        % paired t-test
        [~, pp_bias(tt,:)] = ttest(bias0_temp,bias1_temp);
        
        
        thresh0_temp = cellfun(@(x) x(:,tt), temporal_thresh0);
        thresh1_temp = cellfun(@(x) x(:,tt), temporal_thresh1);
        % Delta threshold
        delta_thresh_temporal = thresh1_temp - thresh0_temp;
        delta_thresh_temporal_mean(tt,:) = nanmean(delta_thresh_temporal);
        delta_thresh_temporal_sem(tt,:) = nanstd(delta_thresh_temporal) ./ sqrt(sum(~isnan(delta_thresh_temporal)));
        % For ANOVA
        delta_thresh4anova(:,:,tt) = delta_thresh_temporal;

        % paired t-test
        [~, pp_thresh(tt,:)] = ttest(thresh0_temp,thresh1_temp);
        
        
        correctrate0_temp = cellfun(@(x) nanmean(x(:,tt)),temporal_correctrate0);      % Average all headings
        correctrate1_temp = cellfun(@(x) nanmean(x(:,tt)),temporal_correctrate1);      % Average all headings
        % Delta Correct rate
        delta_correct_temporal = correctrate1_temp - correctrate0_temp;
        delta_correct_temporal_mean(tt,:) = nanmean(delta_correct_temporal);
        delta_correct_temporal_sem(tt,:) = nanstd(delta_correct_temporal) ./  sqrt(sum(~isnan(delta_correct_temporal)));
        % For ANOVA
        delta_correct4anova(:,:,tt) = delta_correct_temporal;
        
        % paired t-test
        [~, pp_correctrate(tt,:)] = ttest(correctrate0_temp,correctrate1_temp);
        
end

% -- Plotting
set(figure(427),'pos', [50 60 1000, 500], 'Name', 'Delta Bias Dynamic'); clf; hold on;
num_wimdow = size(delta_bias_temporal_mean,1);
xx = 1:num_wimdow;
errorbar(repmat(xx', 1,3), delta_bias_temporal_mean, delta_bias_temporal_sem, 'LineWidth',2);

xlim([0 num_wimdow]); ylims=ylim; 
plot(xlim, [0 0], 'k--', 'linew',1);

sig_bias = pp_bias<=0.05;
for c = 1:3
    if any(sig_bias(:,c))
        text(xx(sig_bias(:,c)), repmat(ylims(2)+(c-1)*0.01, 1,sum(sig_bias(:,c))),...
            '\ast', 'color', con_color(c,:), 'FontWeight','bold');
    end
end

xlabel('Temporal period'); ylabel('\Delta PSE (o)');
SetFigure();


set(figure(428),'pos', [50 60 1000, 500], 'Name', 'Delta Threshold Dynamic'); clf; hold on;
num_wimdow = size(delta_thresh_temporal_mean,1);
xx = 1:num_wimdow;
errorbar(repmat(xx', 1,3), delta_thresh_temporal_mean, delta_thresh_temporal_sem, 'LineWidth',2);

xlim([0 num_wimdow]); ylims=ylim; 
plot(xlim, [0 0], 'k--', 'linew',1);

sig_thresh = pp_thresh<=0.05;
for c = 1:3
    if any(sig_thresh(:,c))
        text(xx(sig_thresh(:,c)), repmat(ylims(2)+(c-1)*0.01, 1,sum(sig_thresh(:,c))),...
            '\ast', 'color', con_color(c,:), 'FontWeight','bold');
    end
end
xlabel('Temporal period'); ylabel('\Delta Threshold (o)');
SetFigure();


set(figure(429),'pos', [50 60 1000, 500], 'Name', 'Delta Correct Rate Dynamic'); clf; hold on;
num_wimdow = size(delta_correct_temporal_mean,1);
xx = 1:num_wimdow;
errorbar(repmat(xx', 1,3), delta_correct_temporal_mean, delta_correct_temporal_sem, 'LineWidth',2);

xlim([0 num_wimdow]); ylims=ylim; 
plot(xlim, [0 0], 'k--', 'linew',1);

sig_correct = pp_correctrate<=0.05;
for c = 1:3
    if any(sig_correct(:,c))
        text(xx(sig_correct(:,c)), repmat(ylims(2)+(c-1)*0.001, 1,sum(sig_correct(:,c))),...
            '\ast', 'color', con_color(c,:), 'FontWeight','bold');
    end
end
xlabel('Temporal period'); ylabel('\Delta Correct rate');
SetFigure();

% %% Figure saving
% saveas(815, 'Population_bar', 'png');    saveas(815, 'Population_bar', 'svg');
% saveas(817, 'Session_comparison_scatter', 'png');    saveas(817, 'Session_comparison_scatter', 'svg');
% saveas(818, 'Condition_comparison_scatter', 'png');    saveas(818, 'Condition_comparison_scatter', 'svg');
% saveas(1018, 'Bias&Thresh_Comparison_bar', 'png');    saveas(818, 'Bias&Thresh_Comparison_bar', 'svg');
% saveas(1116, 'Across-condition Normalized dPSE Comparison', 'png');    saveas(1116, 'Across-condition Normalized dPSE Comparison', 'svg');
% saveas(31, 'dPSE Comparison Across Conditions', 'png');  saveas(31, 'dPSE Comparison Across Conditions', 'svg');
% 
end
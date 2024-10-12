% =========== Session-by-Session Comparison =======
% In this script, comparing inactivation sessions
% before and after chemical manipulation

function InactiveCompare(~)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======= Change the file directory===========

%%%%%%% == Inactivation Sessions == %%%%%%%%%%%%%%%
% directory 1: pre_inactivation control sessions
% directory 2: inactivation sessions

direct_path{1}= 'D:\Paper_rawdata\Raw data\CN\Chemical inactivation\m13_m15_chemical\CN\Mus_pre';
direct_path{2}= 'D:\Paper_rawdata\Raw data\CN\Chemical inactivation\m13_m15_chemical\CN\Mus_0h';

direct_path{1}= 'D:\Paper_rawdata\Raw data\CN\Chemical inactivation\m13_m15_chemical\CN\Sch_pre';
direct_path{2}= 'D:\Paper_rawdata\Raw data\CN\Chemical inactivation\m13_m15_chemical\CN\Sch_post';

direct_path{1}= 'D:\Paper_rawdata\Raw data\CN\Chemical inactivation\m13_m15_chemical\CN\saline_pre';
direct_path{2}= 'D:\Paper_rawdata\Raw data\CN\Chemical inactivation\m13_m15_chemical\CN\saline_post';


direct_path{1}= 'D:\Paper_rawdata\Raw data\CN\Chemical inactivation\m13_m15_chemical\FEF_LIP\LIP\Mus_pre';
direct_path{2}= 'D:\Paper_rawdata\Raw data\CN\Chemical inactivation\m13_m15_chemical\FEF_LIP\LIP\Mus_post';

direct_path{1}= 'D:\Paper_rawdata\Raw data\CN\Chemical inactivation\m13_m15_chemical\FEF_LIP\FEF\Mus_pre';
direct_path{2}= 'D:\Paper_rawdata\Raw data\CN\Chemical inactivation\m13_m15_chemical\FEF_LIP\FEF\Mus_post';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Which directory is experiment group (0h after injection)
experiment_group = 2;
%% Pre-determined
% PSE and Threshold are two variables to be compared
variab{1} = 'Bias_psy';
variab{2} = 'Thresh_psy';
variab{3} = 'Correct_rate';    % Also cate correct rate in inactivation task


set(0, 'defaultFigureRenderer', 'painters')
set(0, 'defaultFigureRendererMode', 'Manual')
set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46]/255);

condi_color = mat2cell([41 89 204; 248 28 83; 14 153 46]/255,ones(3,1));
condi_linestyle = {'b-', 'r-', 'g-'};


%% Extract Data
for di =1: length(direct_path)    % index of directory
    cd(direct_path{di});
    
    % FEF and LIP Data of Fara is from right hemisphere
    if_FEF_LIP = ~isempty(strfind(direct_path{1}, 'FEF_LIP'));
    
    fileName_list{di} = dir('*_Psycho.mat');
    
    for sess = 1:length(fileName_list{di})
        
        load(fileName_list{di}(sess).name);
        if result.repetitionN < 5
            continue;
        end
        
        file_name{di,sess} = result.FILE;
        if_Fara = str2num(result.FILE(2:3)) == 13;
        
        PSE{di,sess} = cell2mat(result.Bias_psy);
        if if_FEF_LIP && if_Fara
            PSE{di,sess} = - PSE{di,sess};
        end
        thresh{di,sess} = cell2mat(result.Thresh_psy);
        Correct_Rate{di,sess} = cell2mat(result.correct_rate);
        raw_data{di,sess} = result.raw;     % For psychometric function fitting and statistical test
        
        
        % Added by ZZ @ 20230722
        % Temporal dynamics of behavior change within the block
        % to see the adaptation
        bias_shift{di,sess} = reshape(cell2mat(result.psy_bias_shift),[], 3);
        rightward_prop_shift{di,sess} = result.right_choice_prop_shift';
        
        % To confirm that the inactivation does not change monkeys' strategy
        beta_with_hist{di,sess} = cell2mat(result.beta_with_hist);     % coefficients of logistic regression: Bias, Slope, Rewarded previous trial, Not rewarded previous trial
        p_value_with_hist{di,sess} = cell2mat(result.p_value_with_hist);
        prop_win_switch{di,sess} = result.prob_cur_cond_win_switch;
        prop_should_win_switch{di,sess} = result.prob_should_win_switch;
        prop_lose_stay{di,sess} = result.prob_cur_cond_lose_stay;
        prop_should_lose_stay{di,sess} = result.prob_should_lose_stay;
    end
end

%%  Align data from different directories
% Which directory is experiment group
% [~, experiment_group] = min(cellfun(@length, fileName_list));
monkey_ind = cellfun(@(x) strfind(x,'m'), file_name, 'UniformOutput', false);
cell_ind = cellfun(@(x) strfind(x,'c'), file_name, 'UniformOutput', false);
run_ind = cellfun(@(x) strfind(x,'r'), file_name, 'UniformOutput', false);

same_sess_ind = repmat(1:size(PSE,2),size(PSE,1),1); % Initiation
for se = 1:length(fileName_list{experiment_group})   % sessions from experiment-session directory
    cell_num = file_name{experiment_group,se}(cell_ind{experiment_group,se}:run_ind{experiment_group,se}-1);
    
    for od = setdiff(1:length(fileName_list),experiment_group)
        % Yong wants to see drug effects after 24h
        %         if monkey_ind == 13 && od == 3
        %             same_cell_ind = cellfun(@(x) strfind(x, num2str(str2num(cell_num)+1)), file_name(od, :), 'UniformOutput', false);  % Pre of the second day is the post of today
        %         else
        same_cell_ind = cellfun(@(x) strfind(x, cell_num), file_name(od, :), 'UniformOutput', false);  % Whether the other directory has corresponding session
        %         end
        
        try
            same_sess_ind(od,se) = find(cellfun(@(x) ~isempty(x), same_cell_ind));  % The index of corresponding session in the other directory
        catch
            same_sess_ind(od,se) = NaN;
            continue;
        end
        
    end
end


%% Rearrange data for plotting
for dl = 1:length(fileName_list)  % directory number
    
    % Initiation for the different repetitions number in different blocks
    arrange_PSE_shift{dl} = nan(length(fileName_list{experiment_group}), max(cellfun(@(x) size(x,1), bias_shift(1,:))),size(bias_shift{1},2));
    arrange_rightward_prop_shift{dl} = nan(length(fileName_list{experiment_group}), max(cellfun(@(x) size(x,1), rightward_prop_shift(1,:))),size(rightward_prop_shift{1},2));
    for se = 1:length(fileName_list{experiment_group})
        
        if isnan(same_sess_ind(dl,se))
            arrange_PSE{dl}(se,:) = nan(1,3);
            arrange_Thresh{dl}(se,:) = nan(1,3);
            arrange_CorrectR{dl}(se,:) =nan(1,3);
            monkey_num(dl,se) = nan(1,1);
        else
            arrange_PSE{dl}(se,:) = PSE{dl,same_sess_ind(dl,se)};
            arrange_Thresh{dl}(se,:) = thresh{dl,same_sess_ind(dl,se)};
            arrange_CorrectR{dl}(se,:) = Correct_Rate{dl,same_sess_ind(dl,se)};
            
            arrange_PSE_shift{dl}(se,1:size(bias_shift{dl,same_sess_ind(dl,se)},1),:) = bias_shift{dl,same_sess_ind(dl,se)};  % sess_num * shift_num * stim_type
            arrange_rightward_prop_shift{dl}(se,1:size(rightward_prop_shift{dl,same_sess_ind(dl,se)},1),:) = rightward_prop_shift{dl,same_sess_ind(dl,se)};
            arrange_beta_with_hist{dl}(se,:,:) = beta_with_hist{dl,same_sess_ind(dl,se)};
            arrange_p_value_with_hist{dl}(se,:,:) = p_value_with_hist{dl,same_sess_ind(dl,se)};
            arrange_prop_win_switch{dl}(se,:) = prop_win_switch{dl,same_sess_ind(dl,se)}';
            arrange_prop_should_win_switch{dl}(se,:) = prop_should_win_switch{dl,same_sess_ind(dl,se)}';
            arrange_prop_lose_stay{dl}(se,:) = prop_lose_stay{dl,same_sess_ind(dl,se)}';
            arrange_prop_should_lose_stay{dl}(se,:) = prop_should_lose_stay{dl,same_sess_ind(dl,se)}';
            
            monkey_num(dl,se) = str2num(file_name{dl,same_sess_ind(dl,se)}(2:3)); % the number of all my monkeys are double-digit
            
        end
        
    end
end

%% Plotting
monkey_name = {'D', 'F'}; monkeyN = [15; 13];
monkey_symbol = {'o', '^'};
monkey_line = {'-', '--'};

%% 1. PSE, Threshold, Correct rate Comparison before and after inactivation
set(figure(501), 'Position', [10 50 1500 400], 'Name', 'Inactivation_Comparison'); clf;
% hs = tight_subplot(3,length(variab),[.1 .1],[.05 .1]);

for st = 1:3  % stim type
    subplot(1,3,st); hold on; 
%     set(gcf, 'CurrentAxes', hs((st-1)*3+1)); hold on;
    
%     [~, pp_2m] = ttest2(arrange_PSE{1}(:,st), arrange_PSE{2}(:,st));
   pp_2m = signrank(arrange_PSE{1}(:,st), arrange_PSE{2}(:,st));

    
    for m = unique(monkey_num(experiment_group,:))
        mi = monkey_num(experiment_group,:) == m;  % monkey_ind
        which_m = find(monkeyN == m);
        
        pse_plot = cell2mat(cellfun(@(x) x(mi,st), arrange_PSE, 'uniform',0));
        plot(1+(which_m-1.5)/5:dl+(which_m-1.5)/5, pse_plot, monkey_line{which_m}, 'color',condi_color{st}, 'linew',2); hold on;
        
        plot(1+(which_m-1.5)/5:dl+(which_m-1.5)/5, pse_plot, monkey_symbol{which_m},...
            'MarkerEdgeColor','k', 'MarkerFaceColor', condi_color{st}, 'MarkerSize',11,'linew',2);
        
        % Paired t-test
%         [~, pp] = ttest2(pse_plot(:,1), pse_plot(:,2));
        pp = signrank(pse_plot(:,1), pse_plot(:,2));
        text(1.1,150+20*(which_m-1), [monkey_name{which_m} ' : ' num2str(pp) ' n = ' num2str(size(pse_plot,1))], 'Color',condi_color{st},  'Units','pixels')
        
    end
    text(1.1,0, ['2Monkeys' ' : ' num2str(pp_2m) ' n = ' num2str(size(arrange_PSE{1},1))], 'Color',condi_color{st},  'Units','pixels')
    
    xlim([1-0.2 dl+0.2])
    xticks(1:1:dl)
    xticklabels({'pre', '0h', 'post'})
    ylabel('PSE')
    
end
SetFigure();

%== Threshold and correct rate 
fig = gcf; fig = fig.Number+1;
set(figure(fig), 'Position', [10 50 1500 900], 'Name', 'Inactivation_Comparison'); clf; hold on;
for st = 1:3  % stim type
    subplot(2,3,st); hold on; 
%     set(gcf, 'CurrentAxes', hs((st-1)*3+2)); hold on;
    
%     [~, pp_2m] = ttest2(arrange_Thresh{1}(:,st), arrange_Thresh{2}(:,st));
    pp_2m = signrank(arrange_Thresh{1}(:,st), arrange_Thresh{2}(:,st));
    for m = unique(monkey_num(experiment_group,:))
        mi = monkey_num(experiment_group,:) == m;  % monkey_ind
        which_m = find(monkeyN == m);
        
        thresh_plot = cell2mat(cellfun(@(x) x(mi,st), arrange_Thresh, 'uniform',0));
        plot(1+(which_m-1.5)/5:dl+(which_m-1.5)/5, thresh_plot, monkey_line{which_m}, 'color',condi_color{st}, 'linew',2); hold on;
        
        plot(1+(which_m-1.5)/5:dl+(which_m-1.5)/5, thresh_plot, monkey_symbol{which_m},...
            'MarkerEdgeColor','k', 'MarkerFaceColor', condi_color{st}, 'MarkerSize',11,'linew',2);
        
        % Paired t-test
%         [~, pp] = ttest2(thresh_plot(:,1), thresh_plot(:,2));
        pp = signrank(thresh_plot(:,1), thresh_plot(:,2));
        text(1.1,150+20*(which_m-1), [monkey_name{which_m} ' : ' num2str(pp) ' n = ' num2str(size(thresh_plot,1))], 'Color',condi_color{st},  'Units','pixels')
    end
    text(1.1,0, ['2Monkeys' ' : ' num2str(pp_2m) ' n = ' num2str(size(arrange_Thresh{1},1))], 'Color',condi_color{st},  'Units','pixels')
    xlim([1-0.2 dl+0.2])
    xticks(1:1:dl)
    xticklabels({'pre', '0h', 'post'})
    ylabel('Threshold')
end


for st = 1:3  % stim type
        subplot(2,3,3+st); hold on; 
    
%     set(gcf, 'CurrentAxes', hs((st-1)*3+3)); hold on;
%     [~, pp_2m] = ttest2(arrange_CorrectR{1}(:,st), arrange_CorrectR{2}(:,st));
    pp_2m = signrank(arrange_CorrectR{1}(:,st), arrange_CorrectR{2}(:,st));

    for m = unique(monkey_num(experiment_group,:))
        mi = monkey_num(experiment_group,:) == m;  % monkey_ind
        which_m = find(monkeyN == m);
        
        correctR_plot = cell2mat(cellfun(@(x) x(mi,st), arrange_CorrectR, 'uniform',0));
        plot(1+(which_m-1.5)/5:dl+(which_m-1.5)/5, correctR_plot, monkey_line{which_m}, 'color',condi_color{st}, 'linew',2); hold on;
        
        plot(1+(which_m-1.5)/5:dl+(which_m-1.5)/5, correctR_plot, monkey_symbol{which_m},...
            'MarkerEdgeColor','k', 'MarkerFaceColor', condi_color{st}, 'MarkerSize',11,'linew',2);
        
        % Paired t-test
%         [~, pp] = ttest2(correctR_plot(:,1), correctR_plot(:,2));
        pp = signrank(correctR_plot(:,1), correctR_plot(:,2));
        text(1.1,150+20*(which_m-1), [monkey_name{which_m} ' : ' num2str(pp) ' n = ' num2str(size(correctR_plot,1))], 'Color',condi_color{st},  'Units','pixels')
    end
    text(1.1,0, ['2Monkeys' ' : ' num2str(pp_2m) ' n = ' num2str(size(arrange_CorrectR{1},1))], 'Color',condi_color{st},  'Units','pixels')
    xlim([1-0.2 dl+0.2])
    xticks(1:1:dl)
    xticklabels({'pre', '0h', 'post'})
    ylabel('Correct Rate')
    
end
SetFigure();


%% Scatter plots
% Added by ZZ @ 20240221
fig = fig +1;
set(figure(fig), 'Position', [10 50 600 600], 'Name', 'Inactivation_Comparison'); clf; hold on;

for st = 1:3  % stim type
    
%     [~, pp_2m] = ttest2(arrange_PSE{1}(:,st), arrange_PSE{2}(:,st));
    pp_2m = signrank(arrange_PSE{1}(:,st), arrange_PSE{2}(:,st));
    
    for m = unique(monkey_num(experiment_group,:))
        mi = monkey_num(experiment_group,:) == m;  % monkey_ind
        which_m = find(monkeyN == m);
        
        pse_plot = cell2mat(cellfun(@(x) x(mi,st), arrange_PSE, 'uniform',0));
        plot(pse_plot(:,1), pse_plot(:,2), monkey_symbol{which_m}, 'color','k', 'MarkerFaceColor',condi_color{st}, 'MarkerSize',11,'LineWidth',2);
        
        % Paired t-test
%         [~, pp] = ttest2(pse_plot(:,1), pse_plot(:,2));
        pp = signrank(pse_plot(:,1), pse_plot(:,2));
        text(1.1,150+20*((which_m-1.5)*3+st), [monkey_name{which_m} ' : ' num2str(pp) ' n = ' num2str(size(pse_plot,1))], 'Color',condi_color{st},  'Units','pixels')
        
    end
    text(1.1,20+10*(st-1), ['2Monkeys' ' : ' num2str(pp_2m) ' n = ' num2str(size(arrange_PSE{1},1))], 'Color',condi_color{st},  'Units','pixels')
    
end
axis square;

% Line of unity
min_scale = min([arrange_PSE{1}(:) ; arrange_PSE{2}(:)]);
max_scale = max([arrange_PSE{1}(:) ; arrange_PSE{2}(:)]);
axis([min_scale-0.1 max_scale+0.1 min_scale-0.1 max_scale+0.1]);
plot([min_scale-0.2 max_scale*1.05], [min_scale-0.2 max_scale*1.05], 'k--');
plot([min_scale-0.2 max_scale*1.05], [0 0], 'k--');
plot([0 0],[min_scale-0.2 max_scale*1.05], 'k--');

xlabel('Pre')
ylabel('Post')
SetFigure();

%% 2. Temporal dynamics of
color = jet;
session_color = round(linspace(1, 64, size(arrange_PSE_shift{1},1)));
session_color = color(session_color', :);

fig=fig+1;
set(figure(fig), 'Position', [10 50 500*length(variab) 900], 'Name', 'PSE temporal changes'); clf;
for st = 1:3
    for dl = 1:2 %length(arrange_PSE_shift)
        subplot(2,3,st+(dl-1)*3); hold on;
        mean_change = nanmean(arrange_PSE_shift{dl}(:,:,st));
        for sess = 1:size(arrange_PSE_shift{dl},1)
            plot(1:size(arrange_PSE_shift{dl},2), arrange_PSE_shift{dl}(sess,:,st),'-', 'color',session_color(sess,:), 'linew',2);
        end
        plot(1:size(arrange_PSE_shift{dl},2), mean_change, 'k--', 'linew',2);
        axis tight
        
        if st ==1
            if dl == 1
                ylabel('No Injection PSE')
            elseif dl ==2
                ylabel('Injection PSE')
            end
        end
        
        if st == 2 && dl==2
            xlabel('Steps')
        end
        
    end
end
SetFigure();

fig=fig+1;
set(figure(fig), 'Position', [10 50 500*length(variab) 900], 'Name', 'Contralateral choice proportion temporal changes'); clf;
for st = 1:3
    for dl = 1: 2 %length(arrange_rightward_prop_shift)
        subplot(2,3,st+(dl-1)*3); hold on;
        mean_change = nanmean(arrange_rightward_prop_shift{dl}(:,:,st));
        for sess = 1:size(arrange_rightward_prop_shift{dl},1)
            plot(1:size(arrange_rightward_prop_shift{dl},2), arrange_rightward_prop_shift{dl}(sess,:,st),'-', 'color',session_color(sess,:), 'linew',2);
        end
        plot(1:size(arrange_rightward_prop_shift{dl},2), mean_change, 'k--', 'linew',2);
        axis tight
        
        if st ==1
            if dl == 1
                ylabel('No Injection Contra. Proportion')
            elseif dl ==2
                ylabel('Injection Contra. Proportion')
            end
        end
        
        if st == 2 && dl==2
            xlabel('Steps')
        end
        
    end
end
SetFigure();


%% 3. Strategy changes ?
% 1). ===== Logistic regression for the influence of reward history of previous
% trial on current-trial choices
fig=fig+1;
set(figure(fig), 'Position', [10 50 1500 900], 'Name', 'Beta of previous trial history'); clf;

dl = 2; %length(arrange_beta_with_hist);
% Reward history of previous trial
for st = 1:3  % stim type
    subplot(2,3,st); hold on;
    
    reward_history = cell2mat(cellfun(@(x) x(:,3,st), arrange_beta_with_hist(1:dl), 'uniform',0)); % Reward history of previous trial
    reward_history_p_value = cell2mat(cellfun(@(x) x(:,3,st), arrange_p_value_with_hist(1:dl), 'uniform',0));   % Whether this term is significant
    sig_ind = reward_history_p_value < 0.05;
    
    for m = unique(monkey_num(experiment_group,:))
        mi = monkey_num(experiment_group,:) == m;
        which_m = find(monkeyN == m);
        
        plot(1+(which_m-1.5)/5 : dl+(which_m-1.5)/5, reward_history(mi',:),...
            monkey_line{which_m}, 'color',condi_color{st}, 'linew',2); hold on;
        
        plot(1+(which_m-1.5)/5 : dl+(which_m-1.5)/5, reward_history(mi',:),...
            monkey_symbol{which_m},  'MarkerEdgeColor',condi_color{st}, 'MarkerFaceColor','none','MarkerSize',11, 'linew',2); hold on;
        
        pp = signrank(reward_history(mi',1), reward_history(mi',2));
        text(1.1,150+20*(which_m-1), [monkey_name{which_m} ' : ' num2str(pp) ' n = ' num2str(size(reward_history(mi',:),1))], 'Color',condi_color{st},  'Units','pixels')

        for dn = 1:dl
            if sum(sig_ind(mi',dn)) > 0
                plot(dn+(which_m-1.5)/5, reward_history(sig_ind(:,dn) & mi',dn), monkey_symbol{which_m},...
                    'MarkerEdgeColor',condi_color{st} , 'MarkerFaceColor',condi_color{st}, 'MarkerSize',11,'linew',2);
            end
        end
    end
    
    % Paired t-test
%     [~, pp] = ttest2(reward_history(:,1), reward_history(:,2));
    pp = signrank(reward_history(:,1), reward_history(:,2));
    text(0.1,10, num2str(pp), 'Color',condi_color{st},  'Units','pixels')
    
    xlim([1-0.2 dl+0.2])
    xticks(1:1:dl)
    xticklabels({'pre', '0h', 'post'})
    ylabel('Beta of Previous Reward')
    
end

% No reward history of previous trial
for st = 1:3  % stim type
    subplot(2,3,3+st); hold on;
    
    reward_history = cell2mat(cellfun(@(x) x(:,4,st), arrange_beta_with_hist(1:dl), 'uniform',0));
    reward_history_p_value = cell2mat(cellfun(@(x) x(:,4,st), arrange_p_value_with_hist(1:dl), 'uniform',0));   % Whether this term is significant
    sig_ind = reward_history_p_value < 0.05;
    
    for m = unique(monkey_num(experiment_group,:))
        mi = monkey_num(experiment_group,:) == m;
        which_m = find(monkeyN == m);
        
        plot(1+(which_m-1.5)/5 : dl+(which_m-1.5)/5, reward_history(mi',:),...
            monkey_line{which_m}, 'color',condi_color{st}, 'linew',2); hold on;
        
        plot(1+(which_m-1.5)/5 : dl+(which_m-1.5)/5, reward_history(mi',:),...
            monkey_symbol{which_m},  'MarkerEdgeColor',condi_color{st}, 'MarkerFaceColor','none','MarkerSize',11, 'linew',2);
        
        
        %         plot(1+(st-1)*0.2:dl+(st-1)*0.2, reward_history, '-', 'color',condi_color{st}, 'linew',2);
        %         plot(1+(st-1)*0.2:dl+(st-1)*0.2, reward_history, 'o',...
        %             'MarkerEdgeColor',condi_color{st} , 'MarkerFaceColor','none' , 'MarkerSize',11,'linew',2);
        %         %     plot(1+(st-1)*0.2:dl+(st-1)*0.2, reward_history, 'o',...
        %         %         'MarkerEdgeColor','k', 'MarkerFaceColor', condi_color{st}, 'MarkerSize',11,'linew',2);
        
        for dn = 1:dl
            if sum(sig_ind(mi',dn)) > 0
                plot(dl+(which_m-1.5)/5, reward_history(sig_ind(:,dn)& mi',dn),monkey_symbol{which_m},...
                    'MarkerEdgeColor',condi_color{st} , 'MarkerFaceColor',condi_color{st}, 'MarkerSize',11,'linew',2);
            end
        end
    end
    
    % Paired t-test
%     [~, pp] = ttest2(reward_history(:,1), reward_history(:,2));
    pp = signrank(reward_history(:,1), reward_history(:,2));
    text(0.1,10, num2str(pp), 'Color',condi_color{st},  'Units','pixels')
    
    xlim([1-0.2 dl+0.2])
    xticks(1:1:dl)
    xticklabels({'pre', '0h', 'post'})
    ylabel('Beta of Previous No Reward')
end
SetFigure();


% ======  2). Win-switch-lose-stay strategy ======
delta_win_SwitchMinusStay = cellfun(@(x,y) x-y, arrange_prop_win_switch, arrange_prop_should_win_switch, 'UniformOutput',0);  % This means win_switch bias - win_stay bias
delta_lose_StayMinusSwitch = cellfun(@(x,y) x-y, arrange_prop_lose_stay, arrange_prop_should_lose_stay, 'UniformOutput',0);  % This means lose_stay bias - lose_switch bias

fig=fig+1;
set(figure(fig), 'Position', [10 50 1500 900], 'Name', 'Win-Switch-Lose-Stay Strategy'); clf;

dl = 2; %length(delta_win_SwitchMinusStay);
% == Plot Win Switch minus Stay
for st = 1:3  % stim type
    subplot(2,3,st); hold on;
    
    delta_win_SMS = cell2mat(cellfun(@(x) x(:,st), delta_win_SwitchMinusStay(1:dl), 'uniform',0));
    
    for m = unique(monkey_num(experiment_group,:))
        mi = monkey_num(experiment_group,:) == m;
        which_m = find(monkeyN == m);
        
        plot(1+(which_m-1.5)/5 : dl+(which_m-1.5)/5, delta_win_SMS(mi',:),...
            monkey_line{which_m}, 'color',condi_color{st}, 'linew',2);
        
        plot(1+(which_m-1.5)/5 : dl+(which_m-1.5)/5, delta_win_SMS(mi',:),...
            monkey_symbol{which_m},  'MarkerEdgeColor','k', 'MarkerFaceColor',condi_color{st},'MarkerSize',11, 'linew',2);
        
        %     plot(1+(st-1)*0.2:dl+(st-1)*0.2, delta_win_SMS, '-', 'color',condi_color{st}, 'linew',2);
        %
        %     plot(1+(st-1)*0.2:dl+(st-1)*0.2, delta_win_SMS, 'o',...
        %         'MarkerEdgeColor', 'k', 'MarkerFaceColor',condi_color{st}, 'MarkerSize',11,'linew',2);
        
        pp = signrank(delta_win_SMS(mi',1), delta_win_SMS(mi',2));
        text(1.1,150+20*(which_m-1), [monkey_name{which_m} ' : ' num2str(pp) ' n = ' num2str(size(delta_win_SMS(mi',:),1))], 'Color',condi_color{st},  'Units','pixels')

    xlim([1-0.2 dl+0.2])
    xticks(1:1:dl)
        xticklabels({'pre', '0h'})
        ylabel('Win Switch - Stay')
        
    end
    % Paired t-test
%     [~, pp] = ttest2(delta_win_SMS(:,1), delta_win_SMS(:,2));
    pp = signrank(delta_win_SMS(:,1), delta_win_SMS(:,2));
    text(0.1,10, num2str(pp), 'Color',condi_color{st},  'Units','pixels')
    
end


% == Plot Lose Stay -Switch
for st = 1:3  % stim type
    subplot(2,3,st+3); hold on;
    
    delta_lose_SMS = cell2mat(cellfun(@(x) x(:,st), delta_lose_StayMinusSwitch(1:dl), 'uniform',0));
    
    for m = unique(monkey_num(experiment_group,:))
        mi = monkey_num(experiment_group,:) == m;
        which_m = find(monkeyN == m);
        
        plot(1+(which_m-1.5)/5 : dl+(which_m-1.5)/5, delta_lose_SMS(mi',:),...
            monkey_line{which_m}, 'color',condi_color{st}, 'linew',2);
        
        plot(1+(which_m-1.5)/5 : dl+(which_m-1.5)/5, delta_lose_SMS(mi',:),...
            monkey_symbol{which_m},  'MarkerEdgeColor','k', 'MarkerFaceColor',condi_color{st},'MarkerSize',11, 'linew',2);
        
        
    xlim([1-0.2 dl+0.2])
    xticks(1:1:dl)
        xticklabels({'pre', '0h'})
        ylabel('Lose Stay - Switch')
        
    end
    % Paired t-test
%     [~, pp] = ttest2(delta_lose_SMS(:,1), delta_lose_SMS(:,2));
    pp = signrank(delta_lose_SMS(:,1), delta_lose_SMS(:,2));
    text(0.1,10, num2str(pp), 'Color',condi_color{st},  'Units','pixels')
end
SetFigure();


end

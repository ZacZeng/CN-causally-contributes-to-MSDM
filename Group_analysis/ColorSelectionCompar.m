% This script aims to compare performance before and after drug injection
% in Color Selection Task (Monkey just need to choose indicated target)

function ColorSelectionCompar(~)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======= Change the file directory===========

%%%%%%%%%%%%%%%== Inactivation Sessions == %%%%%%%%%%
% directory 1: Pre_inactivation control sessions
% directory 2: Inactivation sessions

direct_path{1} = 'D:\Paper_rawdata\Raw_data\CN\Chemical inactivation\m13_m15_color_selection\CN\CN_m13_15_Pre_Inject_ColorSelection\muscimol';
direct_path{2} = 'D:\Paper_rawdata\Raw_data\CN\Chemical inactivation\m13_m15_color_selection\CN\CN_m13_15_After_Inject_ColorSelection\muscimol';

% direct_path{1} = 'D:\Paper_rawdata\Raw_data\CN\Chemical inactivation\m13_m15_color_selection\CN\CN_m13_15_Pre_Inject_ColorSelection\sch';
% direct_path{2} = 'D:\Paper_rawdata\Raw_data\CN\Chemical inactivation\m13_m15_color_selection\CN\CN_m13_15_After_Inject_ColorSelection\sch';
% 
% direct_path{1} = 'D:\Paper_rawdata\Raw_data\CN\Chemical inactivation\m13_m15_color_selection\CN\CN_m13_15_Pre_Inject_ColorSelection\saline';
% direct_path{2} = 'D:\Paper_rawdata\Raw_data\CN\Chemical inactivation\m13_m15_color_selection\CN\CN_m13_15_After_Inject_ColorSelection\saline';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Which directory is experiment group
experiment_group = 2; 

%% Pre-determined
variab{1} = 'Correct rate';
variab{2} = 'Contralateral rate';
variab{3} = 'RT contral';   % Added by ZZ @20240717
variab{4} = 'RT ipsi'; 

set(0, 'defaultFigureRenderer', 'painters')
set(0, 'defaultFigureRendererMode', 'Manual')
set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46]/255);

condi_color = mat2cell([41 89 204; 248 28 83; 14 153 46]/255,ones(3,1));

%% Extract Data
for di = 1:length(direct_path)
    cd(direct_path{di});
    
    % FEF and LIP Data of Fara is from right hemisphere
    if_FEF_LIP = ~isempty(strfind(direct_path{1}, 'FEF_LIP'));
    
    fileName_list{di} = dir('*.mat'); 
    
    for sess = 1:length(fileName_list{di})
        load(fileName_list{di}(sess).name);
        
        file_name{di,sess} = fileName_list{di}(sess).name(1:end-4); 
        if_Fara = str2num(file_name{di,sess}(2:3)) == 13;

        Correct_Rate{di,sess} = result.correct_rate;
        Contralateral_rate{di,sess} = result.rightward_rate;  % I inject drug into left hemisphere for both monkeys
        
        % Added by ZZ @20240717
        RT_contra_sac{di,sess} = result.RT_right_sac; 
        RT_ipsi_sac{di,sess} = result.RT_left_sac;
        
        if if_FEF_LIP && if_Fara
            Contralateral_rate{di,sess} = 1-Contralateral_rate{di,sess};
            RT_contra_sac{di,sess} = result.RT_left_sac;
            RT_ipsi_sac{di,sess} = result.RT_right_sac; 
        end
    end
end

%% Align data from different directories
cell_ind =cellfun(@(x) strfind(x,'c'), file_name, 'UniformOutput', 0); 
run_ind = cellfun(@(x) strfind(x,'r'), file_name, 'UniformOutput', 0);

% Initiation
same_sess_ind = repmat(1:length(file_name(experiment_group,:)), size(file_name,1),1); 
for se = 1:length(fileName_list{experiment_group})
    cell_num =  file_name{experiment_group,se}(cell_ind{experiment_group,se}+1:run_ind{experiment_group,se}-1);
    
    for od = setdiff(1:length(fileName_list), experiment_group)  % For other directories 
        same_cell_ind = cellfun(@(x) strfind(x, ['c' cell_num]), file_name(od, :), 'UniformOutput', false);  % Whether the other directory has corresponding session
        
        try 
            same_sess_ind(od, se) = find(cellfun(@(x) ~isempty(x), same_cell_ind));   % The index of corresponding session in the other directories
        catch
            same_sess_ind(od,se) = nan; 
            continue;
        end
        
    end
end

%% Rearrange data for plotting
for dl = 1:length(fileName_list)   % Directory num
    for se = 1:length(fileName_list{experiment_group})  % Experiment session
        if isnan(same_sess_ind(dl,se))
            arrange_CorrectR(dl,se) = nan;
            arrange_ContralaterR(dl,se) = nan;
            RT_contra(dl,se) = nan;
            RT_ipsi(dl,se) = nan; 
        else
            arrange_CorrectR(dl,se) = Correct_Rate{dl,se};
            arrange_ContralaterR(dl,se) = Contralateral_rate{dl,se};
            
            RT_contra(dl,se) = median(RT_contra_sac{dl,se});
            RT_ipsi(dl,se) = median(RT_ipsi_sac{dl,se}); 
        end
        
    end
end

%% Plotting
set(figure, 'Position', [10 50 500*length(variab) 400], 'name', 'Color Selection Comparison'); clf;

for v = 1:length(variab)
    subplot(1,length(variab),v); hold on;
    if v == 1   % Correct rate
        
        % Median and 95% CI
        medi = nanmedian(arrange_CorrectR, 2);
        low_ci = prctile(arrange_CorrectR, 2.5, 2);
        high_ci = prctile(arrange_CorrectR, 97.5, 2);

        plot(1:size(arrange_CorrectR,1) , arrange_CorrectR, '-', 'color', [0.5 0.5 0.5], 'linew',2);
        
        plot(1:size(arrange_CorrectR,1) , arrange_CorrectR, 'o', 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize',11,'linew',2); 
        
        plot(1:size(arrange_CorrectR,1), medi, 'o',  'MarkerEdgeColor', 'k',  'MarkerFaceColor', 'k', 'MarkerSize',15,'linew',2)
        plot(1:size(arrange_CorrectR,1), [low_ci(1); low_ci(2)], 'r.', 'MarkerSize',10);
        plot(1:size(arrange_CorrectR,1),  [high_ci(1); high_ci(2)], 'r.', 'MarkerSize',10);
        
        % Wilcoxon sign rank test
        pp = signrank(arrange_CorrectR(1,:), arrange_CorrectR(2,:)); 
                
        ylabel(variab{v});
        ylim([0.96 1])
        yticks(0.96:0.01:1)
        
        text(1, 0.982, num2str(pp))
        
    elseif v==2   % Contralateral rate
        
        % Median and 95% CI
        medi = nanmedian(arrange_ContralaterR, 2);
        low_ci = prctile(arrange_ContralaterR, 2.5, 2);
        high_ci = prctile(arrange_ContralaterR, 97.5, 2);

        plot(1:size(arrange_ContralaterR,1) , arrange_ContralaterR, '-', 'color', [0.5 0.5 0.5], 'linew',2);
        
        plot(1:size(arrange_ContralaterR,1) , arrange_ContralaterR, 'o', 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize',11,'linew',2);
        
        
        plot(1:size(arrange_ContralaterR,1), medi, 'o',  'MarkerEdgeColor', 'k',  'MarkerFaceColor', 'k', 'MarkerSize',15,'linew',2)
        plot(1:size(arrange_ContralaterR,1), [low_ci(1); low_ci(2)], 'r.', 'MarkerSize',10);
        plot(1:size(arrange_ContralaterR,1),  [high_ci(1); high_ci(2)], 'r.', 'MarkerSize',10);

                % Wilcoxon sign rank test
        pp = signrank(arrange_ContralaterR(1,:), arrange_ContralaterR(2,:));
        
        ylabel(variab{v});
        ylim([0.47 0.53])
        yticks(0.47:0.03:0.53)
        
        text(1, 0.472, num2str(pp))

    elseif v==3   % Contralateral rate
        
        % Median and 95% CI
        medi = nanmedian(RT_contra, 2);
        low_ci = prctile(RT_contra, 2.5, 2);
        high_ci = prctile(RT_contra, 97.5, 2);

        plot(1:size(RT_contra,1) , RT_contra, '-', 'color', [0.5 0.5 0.5], 'linew',2);
        
        plot(1:size(RT_contra,1) , RT_contra, 'o', 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize',11,'linew',2);
        
        
        plot(1:size(RT_contra,1), medi, 'o',  'MarkerEdgeColor', 'k',  'MarkerFaceColor', 'k', 'MarkerSize',15,'linew',2)
        plot(1:size(RT_contra,1), [low_ci(1); low_ci(2)], 'r.', 'MarkerSize',10);
        plot(1:size(RT_contra,1),  [high_ci(1); high_ci(2)], 'r.', 'MarkerSize',10);

                % Wilcoxon sign rank test
        pp = signrank(RT_contra(1,:), RT_contra(2,:));
        
        ylabel(variab{v});
%         ylim([0.47 0.53])
%         yticks(0.47:0.03:0.53)
        
        text(1,350, num2str(pp))

    elseif v==4   % Contralateral rate
        
        % Median and 95% CI
        medi = nanmedian(RT_ipsi, 2);
        low_ci = prctile(RT_ipsi, 2.5, 2);
        high_ci = prctile(RT_ipsi, 97.5, 2);

        plot(1:size(RT_ipsi,1) , RT_ipsi, '-', 'color', [0.5 0.5 0.5], 'linew',2);
        
        plot(1:size(RT_ipsi,1) , RT_ipsi, 'o', 'MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerSize',11,'linew',2);
        
        
        plot(1:size(RT_ipsi,1), medi, 'o',  'MarkerEdgeColor', 'k',  'MarkerFaceColor', 'k', 'MarkerSize',15,'linew',2)
        plot(1:size(RT_ipsi,1), [low_ci(1); low_ci(2)], 'r.', 'MarkerSize',10);
        plot(1:size(RT_ipsi,1),  [high_ci(1); high_ci(2)], 'r.', 'MarkerSize',10);

                % Wilcoxon sign rank test
        pp = signrank(RT_ipsi(1,:), RT_ipsi(2,:));
        
        ylabel(variab{v});
%         ylim([0.47 0.53])
%         yticks(0.47:0.03:0.53)
        
        text(1, 350, num2str(pp))

    end
    xticks(1:size(arrange_CorrectR,1)); xticklabels({'pre', '0 h'}); 

end
title(['N = ' num2str(sum(~isnan(arrange_ContralaterR')))]);
SetFigure();
end
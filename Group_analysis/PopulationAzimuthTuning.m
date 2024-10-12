% 1DAT population and cell-cell comparison with choice preference in Heading task

clear all; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======= Change the file directory===========

cd 'D:\Paper_rawdata\Raw data\CN\Recordings\CN_m15_1DAT'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-determined
set(0, 'defaultFigureRenderer', 'painters')
set(0, 'defaultFigureRendererMode', 'Manual')
set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46]/255);

condi_color = mat2cell([41 89 204; 248 28 83; 14 153 46]/255,ones(3,1));

%% Extract Data
fileName_list = dir('*_AzimuthTuning.mat');
counter = 1;
for c = 1:length(fileName_list)   % cell num
    load(fileName_list(c).name);
    
    if result.mean_repetitionN < 5 || length(result.unique_stim_type) < 3
        pref_azi(counter, :) = nan(1,3);
        DDI(counter, :) = nan(1,3);
        HTI(counter, :) = nan(1,3);
        p_anova(counter, :) = nan(1,3);
        Dprime(counter,:) = nan(1,3);

        counter = counter + 1;
        continue;
    end
    
%     if counter == 286
%         keyboard;
%     end
    
    pref_azi(counter, :) = result.az;
    DDI(counter, :) = result.DDI; 
    HTI(counter, :) = result.HTI_;
    p_anova(counter, :) = result.p_1D;
    Dprime(counter,:) = result.Dprime;
    
    counter = counter + 1;
    
end


%% Plotting

% 1. DDI
edges = 0:0.1:1;
bin_center = 0.05:0.1:0.95; 
sig_ind = p_anova<0.05;
mean_DDI = nanmean(DDI); 

set(figure, 'name','Direction disrimination index', 'pos',[10 10 1500 1000]); 
suptitle('DDI');
for c = 1:3  % stim_type 
    N_sig(c,:) = histcounts(DDI(sig_ind(:,c),c), edges); 
    N_nonsig(c,:) = histcounts(DDI(~sig_ind(:,c),c), edges); 
    
    % Normalization (Frequency)
    N_sig_norm(c,:) = N_sig(c,:) / sum(~isnan(DDI(:,c))); 
    N_nonsig_norm(c,:) = N_nonsig(c,:) / sum(~isnan(DDI(:,c))); 
    
    N_count(c,:) = histcounts(DDI(:,c), edges); 
    N_count_norm(c,:) = N_count(c,:) / sum(~isnan(DDI(:,c))); 
    
    subplot(2, 3, c); hold on; 
    bar(bin_center, N_count(c,:)', 'FaceColor', 'none', 'EdgeColor',condi_color{c}, 'LineWidth',1.5); 
    bar(bin_center, N_sig(c,:)', 'FaceColor', condi_color{c}, 'EdgeColor',condi_color{c}, 'LineWidth',1.5); 
    xticks(bin_center(1:3:end))
    
    % Normalization 
    subplot(2, 3, c+3); hold on; 
    bar(bin_center, N_count_norm(c,:)', 'FaceColor', 'none', 'EdgeColor',condi_color{c}, 'LineWidth',1.5); 
    bar(bin_center, N_sig_norm(c,:)', 'FaceColor', condi_color{c}, 'EdgeColor',condi_color{c}, 'LineWidth',1.5); 
    % Mean DDI
    plot([mean_DDI(c) mean_DDI(c)], [0 0.5], '--', 'color',condi_color{c}, 'LineWidth',1.5)
    
    title(['n = ' num2str(sum(~isnan(DDI(:,c)))) ' sig-prop = ' num2str(sum(sig_ind(:,c)/sum(~isnan(DDI(:,c))))) ' mean = ' num2str(mean_DDI(c))])
    
    xticks(bin_center(1:3:end)); yticks(0:0.1:0.5);
    xlabel('DDI')
    
end
SetFigure();

figure; hold on;
plot3(DDI(:,1), DDI(:,2), DDI(:,3), 'bo')
plot3([0 1], [0 1], [0 1], 'k-')

p_anova = anova1(DDI)
[corref, p_corr] = corr(DDI, 'rows', 'complete')


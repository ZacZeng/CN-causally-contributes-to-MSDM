% Construct dataset for dPCA (N*S*M*T*R)
% addapted from dpca_demo.m

% trialNum: N * S * M  (0 heading is excluded, only correct trials)
% firingRates: N * S * M * T * maxTrialNum (20 repetitions)
% firingRatesAverage: N * S * M * T

% N is the number of neurons
% S is the number of headings
% M is the number of conditions (ves, vis & comb)
% T is the number of time-points (note that all the trials should have the
% same length in time !)

% trialNum: number of trials for each neuron in each S, M condition (is
% usually different for different conditions and differernt sessions)

% firingRates: all single-trial data together, massive array. Here
% maxTrialNum is the maximum value in trialNum. E.g. if the number of
% trials per condition varied between 1 and 20, then maxTrialNum = 20. For
% the neurons and conditions with less trials, fill remaining entries in
% firingRates with zeros or nans.

% firingRatesAverage: average of firingRates over trials (5th dimension).
% If the firingRates is filled up with nans, then it's simply
%    firingRatesAverage = nanmean(firingRates, 5);
% if it's filled up with zeros (as is convenient if it's stored on hard
% drive as a sparse matrix), then
%     firingRatesAverage = bsxfun(@s, mean(firingRates,5), size(firingRate, 5)./trialNum);

function dPCA_RNN(group_data, sort_id, stimulus_type)
% sort_id, HeadingDis_cum_PSTH, determine which PSTH extracted
% stimulus_type, range from 1 to 3, only when sort_id == 2, the
% stimulus_type is required

% Pre-definition
ifSimultaneousRecording = 0;

colors = [41 89 204; 248 28 83; 14 153 46]/255;
% heading_colors = [0 0 0; 0.3137 0.1961 0.1249; 0.6274 0.3921 0.2497; 0.9363 0.5851 0.3726; 1.0000 0.7812 0.4975];
stimtype_name = {'Vestibular'; 'Visual'; 'Combined'};

% == Figure default
set(0, 'defaultFigureRenderer', 'painters')
set(0, 'defaultFigureRendererMode', 'Manual')
% if sort_id == ALL_CorrectCHOICE || sort_id == All_Choice
set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46]/255);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======= Change the file directory===========

% load D:\Paper_rawdata\Raw_data\RNN\RNN_NMI_A4s_Vdelay_N0.40.35_H_TPNM_N256n0.1dt0.05_C1dcondition_data_list.mat
% load D:\Paper_rawdata\Raw_data\RNN\RNN_NMI_V4s_Vdelay_N0.40.4_H_TPNM_N256n0.1dt0.05_C1d2condition_data_list.mat
load D:\Paper_rawdata\Raw_data\RNN\RNN_NMI_Vdelay_N0.40.4_H_TPNM_N256n0.1dt0.05_C1dcondition_data_list.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = 256;    % number of neurons
ST = 3;                                  % 3 stim_type (vestibular, visual, combined)
D =2;                                   % number of dicisions


%% Extract PSTH
time_markers{1} = [1 40];
rate_ts = 1:40;

%% Rearrange dataset
% Align choice according to PREF and non-PREF
if exist('grand_auc')
    for N = 1:size(data_list{1}, 3)
        
        PREF_NULL = 1:2;
        if pref_right(N) > 0  % Prefer rightward choice, but the current data structure is left choice first
            PREF_NULL = 2:-1:1;
        end
        
        for cho = 1:2
            for st = 1:ST
                PSTH_raw{N,(cho-1)*3 + st}= squeeze(data_list{st,PREF_NULL(cho)}(:,:,N));  % N*6; 3 Modality (PREF) + 3 Modality (NULL)
            end
        end
    end
    
else
    for N = 1:size(data_list{1}, 3)
        for cho = 1:2
            for st = 1:ST
                PSTH_raw{N,(cho-1)*3 + st}= squeeze(data_list{st,cho}(:,:,N));  % N*6; 3 Modality (PREF) + 3 Modality (NULL)
            end
        end
    end
end

reps = cellfun(@(x) size(x,1), PSTH_raw);

R = max(max(reps));   % Repetition is the largest trial number
T = size(PSTH_raw{1}, 2); % number of time bins

% Append conditions without enough trials with NaN
need_append = R - reps;
nan_append = cellfun(@(x) nan(x,T), mat2cell(need_append, ones(size(need_append,1),1),ones(size(need_append,2),1)), 'UniformOutput',0);
appended_raw_PSTH = cellfun(@(x,y) cat(1,x,y), PSTH_raw, nan_append, 'UniformOutput',0);

% Reorgnize matrix data to N * ST * D* T* R
PSTH = double(permute(reshape(cell2mat(appended_raw_PSTH),R,N,T,ST,D),[2 4 5 3 1]));
trialNum = reshape(reps,N,ST,D);   % Real trial number of each condition

% Important parameters
firingRates = PSTH;

% Averaged PSTH
firingRatesAverage = nanmean(firingRates, 5);

%% Define parameter grouping
% *** Do not change this if you don't know what you are doing! ***
% firingRates array has [ N S C T R] size; here we ignore the 1st dimension
% (neurons), i.e. we have the following parameters:
%   1 - Condition
%   2 - Decision
%   3 - time
% So there are three pairwise interactions:
%  [1 3] - Condition/time interaction
%  [2 3] - decision/time interaction
%  [1 2] - Condition/decision interaction
% And one three-way interaction:
%  [1 2 3] - rest
% I group stimulus with stimulus/time interaction etc.:

% if sort_id == ALL_CorrectCHOICE || sort_id ==All_Choice
combinedParams = {{1,[1 3]}, {2,[2,3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Modality', 'Choice', 'Condition-independent', 'M/C Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

% Time events of interest ( stimulus onset or saccade onset)
% They are marked on the plots with vertical lines
time = rate_ts;    % each neuron has the same time window
timeEvents = time_markers{1};

%% Step 1: PCA of the dataset

X = firingRatesAverage(:,:);
X = bsxfun(@minus, X, mean(X,2));

[W_pca,~,~] = svd(X, 'econ');
W_pca = W_pca(:, 1:20);
% W = W(:, 1:15);

% % minimal plotting
% dpca_plot(firingRatesAverage, W, W, @dpca_plot_default);
% SetFigure(12);

% computing explained variance
explVar_pca = dpca_explainedVariance(firingRatesAverage, W_pca, W_pca, ...
    'combinedParams', combinedParams);

% a bit more informative plotting
dpca_plot(firingRatesAverage, W_pca, W_pca, @dpca_plot_default, ...
    'explainedVar', explVar_pca, ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours);
SetFigure(12);

%% Step 2: PCA in each marginalization separately

dpca_perMarginalization(firingRatesAverage, @dpca_plot_default, ...
    'combinedParams', combinedParams);
SetFigure(12);

%% Step 3: dPCA with regularization

% This function takes some minutes to run. It will save the computations
% in a .mat file with a given name. Once computed, you can simply load
% lambdas out of this file:
%   load('tmp_optimalLambdas.mat', 'optimalLambda')

% Please note that this now includes noise covariance matrix Cnoise which
% tends to provide substantial regularization by itself (even with lambda set
% to zero).
% if sort_id == ALL_CorrectCHOICE
optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
    'combinedParams', combinedParams, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 10);  % increase this number to ~10 for better accuracy

SetFigure(12);

Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
    firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);

[W_dpca_regu,V_dpca_regu,whichMarg_dpca_regu] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);

explVar_dpca_regu = dpca_explainedVariance(firingRatesAverage, W_dpca_regu, V_dpca_regu, ...
    'combinedParams', combinedParams);

dpca_plot(firingRatesAverage, W_dpca_regu, V_dpca_regu, @dpca_plot_default, ...
    'explainedVar', explVar_dpca_regu, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg_dpca_regu,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);
SetFigure(12);


%% PCA and trajectory reconstruction
% @20220418
% Changed to follow regularazation dPCA
% For plotting time markers
plotInt = 8;  % plot intervals, in ms
PCA_time_ind = true(size(rate_ts));  % only care fome stimulus onset to offset
PCA_time = rate_ts(PCA_time_ind);
plotInd = round(linspace(1, length(PCA_time), plotInt));
% We reconstruct 3-d trajectory
denoised_dim = 3;

select_dpcs = W_dpca_regu(:,(whichMarg_dpca_regu==1| whichMarg_dpca_regu==2 |whichMarg_dpca_regu==4));


% [Q,~] = qr(select_dpcs);
% weights_dpca_pca = Q(:,1:3);
Q = pca(select_dpcs);
Q = select_dpcs *  Q(:,1:3);


% Projecting the raw data onto the first several eigenvectors
for i = 1:2  % before or after demixed
    if i == 1
        W = W_pca;
    else
        W= Q;
    end
    
    for dim = 1:denoised_dim
        %         if sort_id == CORRECT_ANGLE && isempty(stimulus_type)
        PCA_projPC_temp = reshape(W(:,dim)' * X, size(firingRatesAverage,2), size(firingRatesAverage,3), size(firingRatesAverage,4), []);
        PCA_projPC{i,dim} = PCA_projPC_temp(:,:,:,:);
        
    end
    % Euclidean distance
    Distance_choice_1d{i} = cellfun(@(x) sqrt(squeeze(x(:,1,:) - x(:,2,:)).^2), PCA_projPC(i,:), 'UniformOutput',0);
    Distance_choice_3d{i} = cell2mat(cellfun(@(x) sqrt(x{1}.^2+x{2}.^2+x{3}.^2), Distance_choice_1d(i), 'UniformOutput',0));
    
    % Modality distance
    con_pair = {[1 2]; [1 3]; [2 3]};
    for cp = 1:length(con_pair)
        Distance_modality_1d{i,cp} = cellfun(@(x) abs(squeeze(x(con_pair{cp}(1),:,:) - x(con_pair{cp}(2),:,:))), PCA_projPC(i,:), 'UniformOutput',0);
        Distance_modality_3d{i,cp} = cell2mat(cellfun(@(x) sqrt(x{1}.^2+x{2}.^2+x{3}.^2), Distance_modality_1d(i,cp), 'UniformOutput',0));
    end
end
% pair_color = [0 0 0; 153 153 153; 166 124 82] / 256;
pair_color = [0 0 0; 0 200 256; 256 0 200] / 256;

% == Plotting
DR = {['PCA']; ['dPCA']};  % dimension reduction method
figN = gcf; figN = figN.Number +1;
for i = 1:2    % before or after demixed
    set(figure(figN),'pos',[18 170 898 786],'name',['3-d Population Dynamics ' DR{i}]); clf; hold on; figN=figN+1;
    
    for st = 1:ST
        % Start point of Choosing PREF
        h_perf(st) = plot3(PCA_projPC{i,1}(st,1,1), PCA_projPC{i,2}(st,1,1),PCA_projPC{i,3}(st,1,1),...
            'o', 'color','k', 'MarkerSize',20, 'MarkerFaceColor','k');
        % Start point of Choosing NULL
        h_null(st) = plot3(PCA_projPC{i,1}(st,2,1), PCA_projPC{i,2}(st,2,1),PCA_projPC{i,3}(st,2,1),...
            'o', 'color','k', 'MarkerSize',20, 'linew',3);
        
        % PREF-choosing trajectory
        plot3(squeeze(PCA_projPC{i,1}(st,1,:)), squeeze(PCA_projPC{i,2}(st,1,:)),squeeze(PCA_projPC{i,3}(st,1,:)),...
            '-', 'color',colors(st,:), 'linew',1.5);
        % Null-choosing trajectory
        plot3(squeeze(PCA_projPC{i,1}(st,2,:)), squeeze(PCA_projPC{i,2}(st,2,:)),squeeze(PCA_projPC{i,3}(st,2,:)),...
            ':', 'color',colors(st,:), 'linew',1.5);
        
        % -- Time markers (intermediate point)
        colorsHsv = repmat(rgb2hsv(colors(st,:)), length(plotInd),1);
        colorsHsv(:,2) = linspace(0.2,1,length(plotInd));
        if st == 3
            colorsHsv(:,3) = linspace(.9,colorsHsv(1,3),length(plotInd));
        end
        colorsRGB = flipud(hsv2rgb(colorsHsv));
        
        for pp = 2:length(plotInd)
            plot3(PCA_projPC{i,1}(st,1,plotInd(pp)), PCA_projPC{i,2}(st,1,plotInd(pp)),PCA_projPC{i,3}(st,1,plotInd(pp)),...
                'o', 'color',colorsRGB(pp,:), 'MarkerFaceColor',colorsRGB(pp,:), 'linew',0.1, 'MarkerSize',15);
            
            plot3(PCA_projPC{i,1}(st,2,plotInd(pp)), PCA_projPC{i,2}(st,2,plotInd(pp)),PCA_projPC{i,3}(st,2,plotInd(pp)),...
                'o', 'color',colorsRGB(pp,:), 'MarkerFaceColor','none', 'linew',2, 'MarkerSize',15);
        end
        xlabel('PC1');   ylabel('PC2');   zlabel('PC3');
        grid on;
        %             axis off;
        
        SetFigure()
        view(10, -15);
    end
    
    
    % == Plot Distance
    % Added @ 20210421
    % -- Distance between choices
    set(figure(figN),'pos',[18 170 1500 300],'name',['Euclidean choice distance ' DR{i}]); clf; hold on; figN=figN+1;
    
    % In 3d subspace
    subplot(1,4,1); hold on;
    for st = 1:ST
        plot(rate_ts(PCA_time_ind), Distance_choice_3d{i}(st,:),  '-', 'color',colors(st,:), 'linew',1.5);
    end
    title('3d subsace'); xlabel('Time (ms)'); ylabel('Distance (a.u.)');
    xlim([min(rate_ts(PCA_time_ind)) max(rate_ts(PCA_time_ind))])
    
    % Each PC
    for pc = 1: length(Distance_choice_1d{i})
        subplot(1,4, pc+1); hold on;
        for st = 1:ST
            plot(rate_ts(PCA_time_ind), Distance_choice_1d{i}{pc}(st,:),  '-', 'color',colors(st,:), 'linew',1.5);
        end
        title(['PC' num2str(pc)]); xlabel('Time (ms)'); ylabel('Distance (a.u.)');
        xlim([min(rate_ts(PCA_time_ind)) max(rate_ts(PCA_time_ind))])
        
    end
    SetFigure();
    legend(stimtype_name, 'location', 'northwest', 'FontSize', 8)
    
    % -- Distance between modalities
    set(figure(figN),'pos',[18 170 1500 300],'name',['Euclidean modality distance ' DR{i}]); clf; hold on; figN=figN+1;
    
    % In 3d subspace
    subplot(1,4,1); hold on;
    for cp = 1:length(con_pair)
        plot(rate_ts(PCA_time_ind), Distance_modality_3d{i,cp}(1,:),  '--', 'color',pair_color(cp,:), 'linew',1.5);  % PREF choice
        plot(rate_ts(PCA_time_ind), Distance_modality_3d{i,cp}(2,:),  ':', 'color',pair_color(cp,:), 'linew',1.5);  % NULL
        plot(rate_ts(PCA_time_ind), (Distance_modality_3d{i,cp}(1,:)+Distance_modality_3d{i,cp}(2,:))/2,  '-', 'color',pair_color(cp,:), 'linew',3);  % Average
    end
    title('3d subsace'); xlabel('Time (ms)'); ylabel('Distance (a.u.)');
    xlim([min(rate_ts(PCA_time_ind)) max(rate_ts(PCA_time_ind))])
    
    % Each PC
    for pc = 1: length(Distance_modality_1d{i})
        subplot(1,4, pc+1); hold on;
        for cp = 1:length(con_pair)
            plot(rate_ts(PCA_time_ind), Distance_modality_1d{i,cp}{pc}(1,:),  '--', 'color',pair_color(cp,:), 'linew',1.5);
            plot(rate_ts(PCA_time_ind), Distance_modality_1d{i,cp}{pc}(2,:),  ':', 'color',pair_color(cp,:), 'linew',1.5);
            plot(rate_ts(PCA_time_ind), (Distance_modality_1d{i,cp}{pc}(1,:)+Distance_modality_1d{i,cp}{pc}(2,:))/2,  '-', 'color',pair_color(cp,:), 'linew',3);
        end
        title(['PC' num2str(pc)]); xlabel('Time (ms)'); ylabel('Distance (a.u.)');
        xlim([min(rate_ts(PCA_time_ind)) max(rate_ts(PCA_time_ind))])
        
    end
    SetFigure();
    legend_txt = {'12 PREF'; '12 NULL';'12 Mean';  '13 PREF'; '13 NULL'; '13 Mean'; '23 PREF'; '23 NULL'; '23 Mean'; };
    legend(legend_txt, 'location', 'northwest', 'FontSize', 6)
    
end

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
%     firingRatesAverage = bsxfun(@times, mean(firingRates,5), size(firingRate, 5)./trialNum);

function dPCA_ZZ(group_data, sort_id, stimulus_type)
% sort_id, HeadingDis_cum_PSTH, determine which PSTH extracted
% stimulus_type, range from 1 to 3, only when sort_id == 2, the
% stimulus_type is required

% Pre-definition
representative_cell = 25;

ifSimultaneousRecording = false;

%%%%%%% Order corresponding to "Sort_Id" in TEMPO_GUI processing %%%%
ALL_CorrectCHOICE = 1; CORRECT_ANGLE = 2; CHOICE_DIFFICULT = 3; OUTCOME = 4; WRONG_ANGLE = 5; CORRECTNWRONG_ANGLE = 6; All_Choice = 7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colors = [41 89 204; 248 28 83; 14 153 46]/255;
heading_colors = [0 0 0; 0.3137 0.1961 0.1249; 0.6274 0.3921 0.2497; 0.9363 0.5851 0.3726; 1.0000 0.7812 0.4975];
stimtype_name = {'Vestibular'; 'Visual'; 'Combined'};

% == Figure default
set(0, 'defaultFigureRenderer', 'painters')
set(0, 'defaultFigureRendererMode', 'Manual')
if sort_id == ALL_CorrectCHOICE || sort_id == All_Choice
    set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46]/255);
    % set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46; 248 28 248; 14 248 248]/255);
elseif ~isempty(stimulus_type) && sort_id == CORRECT_ANGLE || sort_id == CORRECTNWRONG_ANGLE
    set(0,'defaultAxesColorOrder', [0 0 0; 0.3137 0.1961 0.1249; 0.6274 0.3921 0.2497; 0.9363 0.5851 0.3726; 1.0000 0.7812 0.4975]);
    
elseif isempty(stimulus_type) && (sort_id == CORRECT_ANGLE || sort_id == CORRECTNWRONG_ANGLE)
    % Gouki's suggestion: trajectory of three modalities in the same dPCA subspace
    % Colors for different headings in different modalities
    for k = 1:3
        colors_angles{k} = colormap(gray);
        colors_angles{k} = ones(size(colors_angles{k},1),3) - colors_angles{k} .* repmat([1 1 1]-colors(k,:),size(colors_angles{k},1),1);
        colors_angles{k} = colors_angles{k}(round(linspace(20,length(colors_angles{k}),5)),:);
    end
    
    set(0, 'defaultAxesColorOrder', cell2mat(colors_angles'))

end

N = size(group_data,2);    % number of neurons
% T =180;                                       % number of time points
% % H = 4;                                   % 8 headings (excluding 0), 4 pairs
H = 5;                                   % 10 headings (including 0), 5 pairs
ST = 3;                                  % 3 stim_type (vestibular, visual, combined)
% R = 20;                                % maximal number of repetitions
D =2;                                   % number of dicisions
%Because the directions of stimulus and choice are totally related, so
% When only correct trials are included, conditionNum = S * C
% When false trials are included, conditionNum = S * C * D

% if sort_id == CORRECT_ANGLE
%     if  isempty(stimulus_type) || ~ismember(stimulus_type, 1:3)
%         error('Does not offer which stimulus type to be analysed');
%     end
% end


%% Extract PSTH
rate_ts = {group_data(representative_cell).PSTH{1,1,1}.ts, group_data(representative_cell).PSTH{2,1,1}.ts,...
    group_data(representative_cell).PSTH{3,1,1}.ts};

align_markers = group_data(representative_cell).align_markers;
align_offsets_others = group_data(representative_cell).align_offsets_others;
time_markers{1} = [0 mean(align_offsets_others{1})];
time_markers{2} = [0];
time_markers{3} = [0];

unique_heading = unique(group_data(representative_cell).trialInfo(:,2));

PSTH_temp = [group_data(:).PSTH];
PSTH_sort_id = squeeze(PSTH_temp(1,sort_id:7:end,:));
% PSTH = squeeze(PSTH_temp(1,2:7:end,:));     % only want stimulus period, CORRECT_ANGLE trials


%% Rearrange dataset
if sort_id == ALL_CorrectCHOICE || sort_id == All_Choice    % Modality and choice as
    PSTH_sort_id = PSTH_sort_id(cellfun(@(x) ~isempty(x),PSTH_sort_id));
    for c = 1 : ST * D   % number of conditions
        raw_PSTH(:,c) = cellfun(@(x) x.raw{c}, PSTH_sort_id, 'UniformOutput',0);
        reps(:,c) = cellfun(@(x) size(x,1), raw_PSTH(:,c));
    end
    max_reps = max(max(reps));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R = max_reps;
    T =size(PSTH_sort_id{representative_cell}.raw{1},2);                                       % number of time points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    need_append = R - reps;
    nan_append = cellfun(@(x) nan(x,T), mat2cell(need_append, ones(size(need_append,1),1),ones(size(need_append,2),1)), 'UniformOutput',0);
    % Append NaNs to trial conditions without R reps
    appended_raw_PSTH = cellfun(@(x,y) cat(1,x,y), raw_PSTH, nan_append, 'UniformOutput',0);
    
    % Reorgnize matrix data to N * ST * D* T* R
    PSTH = permute(reshape(cell2mat(appended_raw_PSTH),R,N,T,D,ST),[2 5 4 3 1]);
    trialNum = permute(reshape(reps,N,D,ST),[1 3 2]);   % Real trial number of each condition
    
    % Delete sessions with 0 reps in some conditions
    %     [r,~] = find(reps==0);
    % % dpca_optimizeLambda need trialNum of each condition > 1
    % % delete session with only one repetitions in some conditions
    [r,~] = find(reps<=1);
    r = unique(r);
    
    if ~isempty(r)
        PSTH = PSTH(setdiff(1:N,r),:,:,:,:);
        trialNum = trialNum(setdiff(1:N,r),:,:);
    end
    
    % Important parameters
    firingRates = PSTH;
    
    % Averaged PSTH
    firingRatesAverage = nanmean(firingRates, 5);

elseif sort_id == CORRECT_ANGLE   % Sensory and choice in each stim type
    % Drop session without 0 heading
    num_heading = cellfun(@(x) size(x.raw,1), PSTH_sort_id);
    [r,~] = find(num_heading<H*D);
    PSTH_sort_id = PSTH_sort_id(setdiff(1:N,unique(r)),:);
    N = size(PSTH_sort_id,1);
    
    % Drop sessions with more than 20 reps
    % I have nearly no such sessions
    for h = 1 : H*D
        raw_PSTH = cellfun(@(x) x.raw{h}, PSTH_sort_id, 'UniformOutput',0);
        reps(:,:,h) = cellfun(@(x) size(x,1), raw_PSTH);
    end
%     [r,~] = find(reps(:,:)>20);
%     PSTH_sort_id = PSTH_sort_id(setdiff(1:N,unique(r)),:);
%     reps = reps(setdiff(1:N,unique(r)),:,:);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T =size(PSTH_sort_id{representative_cell}.raw{1},2);                                       % number of time points
    N = size(PSTH_sort_id,1);
%     R = 20;
    R = max(max(squeeze(max(reps))));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for h = 1 : H*D
        raw_PSTH = cellfun(@(x) x.raw{h}, PSTH_sort_id, 'UniformOutput',0);
        need_append = R-reps(:,:,h);
        nan_append = cellfun(@(x) nan(x,T), mat2cell(need_append, ones(size(need_append,1),1),ones(size(need_append,2),1)), 'UniformOutput',0);
        
        % Append NaNs to trial conditions without 20 reps
        appended_raw_PSTH(:,:,h) = cellfun(@(x,y) cat(1,x,y), raw_PSTH, nan_append, 'UniformOutput',0);
    end
    
    
    % Data rearrangement
    appended_PSTH = permute(reshape(cell2mat(appended_raw_PSTH), R, N, T, ST, H*D), [2,4,5,3,1]);     % N * H * M * T * maxTrialNum
    
    % Group headings according choosing PREF or NULL
    append_PSTH_pref = appended_PSTH(:,:,1:2:end,:,:);
    reps_pref = reps(:,:,1:2:end);
    append_PSTH_null = appended_PSTH(:,:,2:2:end,:,:);
    reps_null = reps(:,:,2:2:end);
    for st = 1:ST
        PSTH_stimtype{st} = reshape(cat(2,squeeze(append_PSTH_pref(:,st,:,:,:)),squeeze(append_PSTH_null(:,st,:,:,:))), N,H,D,T,R);
        reps_stimtype{st} = reshape(cat(2,squeeze(reps_pref(:,st,:)),squeeze(reps_null(:,st,:))),N,H,D);
    end
    
    % Delete sessions with 0 reps in some conditions
    %     [r,~] = find(reps == 0);
    % % dpca_optimizeLambda need trialNum of each condition > 1
    % % delete session with only one repetitions in some conditions
    [r,~] = find(reps(:,:) <= 1);
    r = unique(r);
    
    if ~isempty(r)
        PSTH_st = cellfun(@(x) x(setdiff(1:N,r),:,:,:,:),PSTH_stimtype, 'UniformOutput',0);
        reps_st = cellfun(@(x) x(setdiff(1:N,r),:,:),reps_stimtype, 'UniformOutput',0);
    end
    
        % We only show 1 stim type, change here to show other stim type
        if ~isempty(stimulus_type)
            st = stimulus_type;
            firingRates = PSTH_st{st};
            trialNum = reps_st{st};
            
            % Averaged PSTH
            firingRatesAverage = nanmean(firingRates, 5);
        end
end
        

if sort_id == CORRECT_ANGLE && isempty(stimulus_type)   % Added by ZZ @ 20230804, suggestions from Gouki
    firingRates = reshape(cell2mat(PSTH_st), size(PSTH_st{1},1), size(PSTH_st{1},2), length(PSTH_st), size(PSTH_st{1},3),...
        size(PSTH_st{1},4), []);
    firingRates = permute(firingRates, [1 2 4 3 5 6]); 
    firingRatesAverage = nanmean(firingRates, 6);
    
    trialNum = reshape(cell2mat(reps_st), size(reps_st{1},1), size(reps_st{1},2), length(reps_st), size(reps_st{1},3)); 
    trialNum = permute(trialNum, [1 2 4 3]); 
end
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

if sort_id == ALL_CorrectCHOICE || sort_id ==All_Choice
    combinedParams = {{1,[1 3]}, {2,[2,3]}, {3}, {[1 2], [1 2 3]}};
    margNames = {'Modality', 'Choice', 'Condition-independent', 'M/C Interaction'};
    margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
    
elseif sort_id == CORRECT_ANGLE && ~isempty(stimulus_type)
    combinedParams = {{1,[1 3]}, {2,[2,3]}, {3}, {[1 2], [1 2 3]}};
    margNames = {'Heading', 'Choice', 'Condition-independent', 'H/C Interaction'};
    margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
    
elseif sort_id == CORRECT_ANGLE && isempty(stimulus_type)
%     combinedParams = {{1,[1 4]}, {2,[2,4]}, {3,[3,4]}, {4}, {[2 3], [2 3 4]}};
%     margNames = {'Heading', 'Modality', 'Choice', 'Condition-independent', 'M/C Interaction'};
%     margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171; 235 19 150]/256;
    
    combinedParams = {{1,[1 4]}, {2,[2,4]}, {3,[3,4]}, {4}};
    margNames = {'Heading', 'Choice', 'Modality', 'Condition-independent'};
    margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;
end

% Time events of interest ( stimulus onset or saccade onset)
% They are marked on the plots with vertical lines
time = rate_ts;    % each neuron has the same time window
timeEvents = time_markers{1};


% %% Check consistency between trialNum and firingRates
% for n = 1:size(firingRates,1)
%     for s = 1:size(firingRates,2)
%         for d = 1:size(firingRates,3)
%             assert(isempty(find(isnan(firingRates(n,s,d,:,1:trialNum(n,s,d))), 1)), 'Something is wrong!')
%         end
%     end
% end


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
    'time', time{1},                        ...
    'timeEvents', timeEvents,               ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours);
SetFigure(12);

%% Step 2: PCA in each marginalization separately

% dpca_perMarginalization(firingRatesAverage, @dpca_plot_default, ...
%     'combinedParams', combinedParams);
% SetFigure(12);

%% Step 3: dPCA without regularization and ignoring noise covariance

% This is the core function.
% W is the decoder, V is the encoder (ordered by explained variance),
% whichMarg is an array that tells you which component comes from which
% marginalization

tic
[W_dpca,V_dpca,whichMarg_dpca] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams);

% [W,V,whichMarg] = dpca(firingRatesAverage, 15, ...
%     'combinedParams', combinedParams);
toc

explVar_dpca = dpca_explainedVariance(firingRatesAverage, W_dpca, V_dpca, ...
    'combinedParams', combinedParams);

dpca_plot(firingRatesAverage, W_dpca, V_dpca, @dpca_plot_default, ...
    'explainedVar', explVar_dpca, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg_dpca,                 ...
    'time', time{1},                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);
SetFigure(12);



%% Step 4: dPCA with regularization

% This function takes some minutes to run. It will save the computations
% in a .mat file with a given name. Once computed, you can simply load
% lambdas out of this file:
%   load('tmp_optimalLambdas.mat', 'optimalLambda')

% Please note that this now includes noise covariance matrix Cnoise which
% tends to provide substantial regularization by itself (even with lambda set
% to zero).
if sort_id == ALL_CorrectCHOICE
    optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
        'combinedParams', combinedParams, ...
        'simultaneous', ifSimultaneousRecording, ...
        'numRep', 10, ...  % increase this number to ~10 for better accuracy
        'filename', 'dPCA_regularization_optimalLambdas_MD_sortid_1.mat');
    
elseif sort_id ==All_Choice
    
    optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
        'combinedParams', combinedParams, ...
        'simultaneous', ifSimultaneousRecording, ...
        'numRep', 10, ...  % increase this number to ~10 for better accuracy
        'filename', 'dPCA_regularization_optimalLambdas_MD_sortid_7.mat');
elseif sort_id == CORRECT_ANGLE && ~isempty(stimulus_type)
    
    optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
        'combinedParams', combinedParams, ...
        'simultaneous', ifSimultaneousRecording, ...
        'numRep', 10, ...  % increase this number to ~10 for better accuracy
        'filename', ['dPCA_regularization_optimalLambdas_SD_' num2str(stimulus_type) '.mat']);
    
elseif sort_id == CORRECT_ANGLE && isempty(stimulus_type)
        optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
        'combinedParams', combinedParams, ...
        'simultaneous', ifSimultaneousRecording, ...
        'numRep', 10, ...  % increase this number to ~10 for better accuracy
        'filename', ['dPCA_regularization_optimalLambdas_SD' '.mat']);

end
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
    'time', time{1},                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);
SetFigure(12);


%% PCA and trajectory reconstruction
% @20220418
% Changed to follow regularazation dPCA 
% For plotting time markers
j=1;
% plotInt = 150;  % plot intervals, in ms
plotInt = 300;  % plot intervals, in ms
PCA_time_ind = rate_ts{j}>= time_markers{j}(1) & rate_ts{j}<= time_markers{j}(2);  % only care fome stimulus onset to offset
PCA_time = rate_ts{j}(PCA_time_ind);
plotPerTimeBin = fix(plotInt / (PCA_time(2) - PCA_time(1)));
plotMinInd = fix(-min(PCA_time)/plotInt)*plotPerTimeBin;      % find the minimal time point (stimulus onset)
plotInd = plotMinInd : plotPerTimeBin : length(PCA_time);
% We reconstruct 3-d trajectory
denoised_dim = 3;

select_dpcs = W_dpca_regu(:,(whichMarg_dpca_regu==1| whichMarg_dpca_regu==2 |whichMarg_dpca_regu==4));
if sort_id == CORRECT_ANGLE && isempty(stimulus_type)
    select_dpcs = W_dpca_regu(:,(whichMarg_dpca_regu==1| whichMarg_dpca_regu==2 |whichMarg_dpca_regu==3));
end

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
        if sort_id == CORRECT_ANGLE && isempty(stimulus_type)
            PCA_projPC_temp = reshape(W(:,dim)' * X, size(firingRatesAverage,2), size(firingRatesAverage,3), size(firingRatesAverage,4), []);
            PCA_projPC{i,dim} = PCA_projPC_temp(:,:,:,PCA_time_ind); 
            
        else
            PCA_projPC_temp = reshape(W(:,dim)' * X, size(firingRatesAverage,2), size(firingRatesAverage,3), []);
            %         PCA_projPC_temp = reshape(W(:,dim)' * X, size(firingRatesAverage,2), size(firingRatesAverage,3), 5, []);
            PCA_projPC{i,dim} = PCA_projPC_temp(:,:,PCA_time_ind);
            %         PCA_projPC{i,dim} = PCA_projPC_temp(:,:,:,PCA_time_ind);
            
        end
    end
    
        if sort_id == ALL_CorrectCHOICE
            % Euclidean distance
            Distance_choice_1d{i} = cellfun(@(x) sqrt(squeeze(x(:,1,:) - x(:,2,:)).^2), PCA_projPC(i,:), 'UniformOutput',0);
            Distance_choice_3d{i} = cell2mat(cellfun(@(x) sqrt(x{1}.^2+x{2}.^2+x{3}.^2), Distance_choice_1d(i), 'UniformOutput',0));
            
            % Modality distance
            con_pair = {[1 2]; [1 3]; [2 3]};
            pair_color =  [0 0 0; 0 200 256; 256 0 200] / 256;
%             pair_color = [0 0 0; 153 153 153; 166 124 82] / 256; 
            for cp = 1:length(con_pair)
                Distance_modality_1d{i,cp} = cellfun(@(x) abs(squeeze(x(con_pair{cp}(1),:,:) - x(con_pair{cp}(2),:,:))), PCA_projPC(i,:), 'UniformOutput',0);
                Distance_modality_3d{i,cp} = cell2mat(cellfun(@(x) sqrt(x{1}.^2+x{2}.^2+x{3}.^2), Distance_modality_1d(i,cp), 'UniformOutput',0));
            end
            
        end
end

DR = {['PCA']; ['dPCA']};  % dimension reduction method
figN = gcf; figN = figN.Number +1;
for i = 1:2    % before or after demixed
    if sort_id == ALL_CorrectCHOICE || sort_id == CORRECT_ANGLE
        % == Plotting
        set(figure(figN),'pos',[18 170 898 786],'name',['3-d Population Dynamics ' DR{i}]); clf; hold on; figN=figN+1;
        
        for st = 1:ST
            
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
            
            % Start point of Choosing PREF
%             h_perf(st) = plot3(PCA_projPC{i,1}(st,1,1), PCA_projPC{i,2}(st,1,1),PCA_projPC{i,3}(st,1,1),...
%                 'o', 'color',colors(st,:), 'MarkerSize',20, 'MarkerFaceColor',colors(st,:));
%             % Start point of Choosing NULL
%             h_null(st) = plot3(PCA_projPC{i,1}(st,2,1), PCA_projPC{i,2}(st,2,1),PCA_projPC{i,3}(st,2,1),...
%                 'o', 'color',colors(st,:), 'MarkerSize',20, 'linew',3);
            h_perf(st) = plot3(PCA_projPC{i,1}(st,1,1), PCA_projPC{i,2}(st,1,1),PCA_projPC{i,3}(st,1,1),...
                'o', 'color','k', 'MarkerSize',20, 'MarkerFaceColor','k');
            % Start point of Choosing NULL
            h_null(st) = plot3(PCA_projPC{i,1}(st,2,1), PCA_projPC{i,2}(st,2,1),PCA_projPC{i,3}(st,2,1),...
                'o', 'color','k', 'MarkerSize',20, 'linew',3);

            
            xlabel('PC1');   ylabel('PC2');   zlabel('PC3'); 
            grid on;
            %             axis off;
            
            SetFigure()
            view(10, -15);
        end
        
        % == Plot Distance 
        % Added @ 20210421
        if sort_id == ALL_CorrectCHOICE
            
                % -- Distance between choices
                set(figure(figN),'pos',[18 170 1500 300],'name',['Euclidean choice distance ' DR{i}]); clf; hold on; figN=figN+1;
                
                % In 3d subspace
                subplot(1,4,1); hold on;  
                for st = 1:ST
                    plot(rate_ts{1}(PCA_time_ind), Distance_choice_3d{i}(st,:),  '-', 'color',colors(st,:), 'linew',1.5);
                end
                title('3d subsace'); xlabel('Time (ms)'); ylabel('Distance (a.u.)');
                xlim([min(rate_ts{1}(PCA_time_ind)) max(rate_ts{1}(PCA_time_ind))])
                
                % Each PC
                for pc = 1: length(Distance_choice_1d{i})  
                    subplot(1,4, pc+1); hold on; 
                    for st = 1:ST
                        plot(rate_ts{1}(PCA_time_ind), Distance_choice_1d{i}{pc}(st,:),  '-', 'color',colors(st,:), 'linew',1.5);
                    end
                    title(['PC' num2str(pc)]); xlabel('Time (ms)'); ylabel('Distance (a.u.)');
                    xlim([min(rate_ts{1}(PCA_time_ind)) max(rate_ts{1}(PCA_time_ind))])
                    
                end
                SetFigure();
                legend(stimtype_name, 'location', 'northwest', 'FontSize', 8)
                
                % -- Distance between modalities
                set(figure(figN),'pos',[18 170 1500 300],'name',['Euclidean modality distance ' DR{i}]); clf; hold on; figN=figN+1;
                
                % In 3d subspace
                subplot(1,4,1); hold on; 
                for cp = 1:length(con_pair)
                    plot(rate_ts{1}(PCA_time_ind), Distance_modality_3d{i,cp}(1,:),  '--', 'color',pair_color(cp,:), 'linew',1.5);  % PREF choice
                    plot(rate_ts{1}(PCA_time_ind), Distance_modality_3d{i,cp}(2,:),  ':', 'color',pair_color(cp,:), 'linew',1.5);  % NULL
                    plot(rate_ts{1}(PCA_time_ind), (Distance_modality_3d{i,cp}(1,:)+Distance_modality_3d{i,cp}(2,:))/2,  '-', 'color',pair_color(cp,:), 'linew',3);  % Average
                end
                title('3d subsace'); xlabel('Time (ms)'); ylabel('Distance (a.u.)');
                xlim([min(rate_ts{1}(PCA_time_ind)) max(rate_ts{1}(PCA_time_ind))])
                
                % Each PC
                for pc = 1: length(Distance_modality_1d{i})
                    subplot(1,4, pc+1); hold on;
                    for cp = 1:length(con_pair)
                        plot(rate_ts{1}(PCA_time_ind), Distance_modality_1d{i,cp}{pc}(1,:),  '--', 'color',pair_color(cp,:), 'linew',1.5);
                        plot(rate_ts{1}(PCA_time_ind), Distance_modality_1d{i,cp}{pc}(2,:),  ':', 'color',pair_color(cp,:), 'linew',1.5);
                        plot(rate_ts{1}(PCA_time_ind), (Distance_modality_1d{i,cp}{pc}(1,:)+Distance_modality_1d{i,cp}{pc}(2,:))/2,  '-', 'color',pair_color(cp,:), 'linew',3);
                    end
                    title(['PC' num2str(pc)]); xlabel('Time (ms)'); ylabel('Distance (a.u.)');
                    xlim([min(rate_ts{1}(PCA_time_ind)) max(rate_ts{1}(PCA_time_ind))])
                    
                end
                SetFigure();
                legend_txt = {'12 PREF'; '12 NULL';'12 Mean';  '13 PREF'; '13 NULL'; '13 Mean'; '23 PREF'; '23 NULL'; '23 Mean'; }; 
                legend(legend_txt, 'location', 'northwest', 'FontSize', 6)
                
                % Normalized comparison, Gu wants it 
                % I computed the AUC
                for cp = 1:length(con_pair)
                    AUC_modality_distance_3d_pref(i,cp) = trapz(rate_ts{1}(PCA_time_ind), Distance_modality_3d{i,cp}(1,:));
                    AUC_modality_distance_3d_null(i,cp) = trapz(rate_ts{1}(PCA_time_ind), Distance_modality_3d{i,cp}(2,:));
                    AUC_modality_distance_3d_mean(i,cp) = trapz(rate_ts{1}(PCA_time_ind), (Distance_modality_3d{i,cp}(1,:)+Distance_modality_3d{i,cp}(2,:))/2);
                end
                
        end
        
    elseif sort_id == 2
        if  isempty(stimulus_type)
            figure('position', [18 170 898 786], 'name', '3-d Population Dynamics (All headings and modalities)'); hold on; 
            
            for st = 1:3
                for hh = 1:H
                    % Start point of Choosing PREF
                    h_perf(hh) = plot3(PCA_projPC{i,1}(hh,1,st,1), PCA_projPC{i,2}(hh,1,st,1),PCA_projPC{i,3}(hh,1,st,1),...
                        'o', 'color',colors_angles{st}(hh,:), 'MarkerSize',20, 'MarkerFaceColor',colors_angles{st}(hh,:));
                    % Start point of Choosing NULL
                    h_null(hh) = plot3(PCA_projPC{i,1}(hh,2,st,1), PCA_projPC{i,2}(hh,2,st,1),PCA_projPC{i,3}(hh,2,st,1),...
                        'o', 'color',colors_angles{st}(hh,:), 'MarkerSize',20, 'linew',3);
                    
                    % PREF-choosing trajectory
                    plot3(squeeze(PCA_projPC{i,1}(hh,1,st,:)), squeeze(PCA_projPC{i,2}(hh,1,st,:)),squeeze(PCA_projPC{i,3}(hh,1,st,:)),...
                        '-', 'color',colors_angles{st}(hh,:), 'linew',1.5);
                    % Null-choosing trajectory
                    plot3(squeeze(PCA_projPC{i,1}(hh,2,st,:)), squeeze(PCA_projPC{i,2}(hh,2,st,:)),squeeze(PCA_projPC{i,3}(hh,2,st,:)),...
                        ':', 'color',colors_angles{st}(hh,:), 'linew',1.5);
                    
                    for pp = 2:length(plotInd)
                        plot3(PCA_projPC{i,1}(hh,1,st,plotInd(pp)), PCA_projPC{i,2}(hh,1,st,plotInd(pp)),PCA_projPC{i,3}(hh,1,st,plotInd(pp)),...
                            'o', 'color',colors_angles{st}(hh,:), 'MarkerFaceColor',colors_angles{st}(hh,:), 'linew',0.1, 'MarkerSize',15);
                        
                        plot3(PCA_projPC{i,1}(hh,2,st,plotInd(pp)), PCA_projPC{i,2}(hh,2,st,plotInd(pp)),PCA_projPC{i,3}(hh,2,st,plotInd(pp)),...
                            'o', 'color',colors_angles{st}(hh,:), 'MarkerFaceColor','none', 'linew',2, 'MarkerSize',15);
                    end
                    
                    grid on;
%                     axis off;

                end
            end
            
        else
        % == Plotting
        figN = gcf; figN = figN.Number +1;
        set(figure(figN),'pos',[18 170 898 786],'name',['3-d Population Dynamics PCA Stimtype = ' num2str(stimulus_type)]); clf; hold on;
        
        for hh = 1:H
            % Start point of Choosing PREF
            h_perf(hh) = plot3(PCA_projPC{i,1}(hh,1,1), PCA_projPC{i,2}(hh,1,1),PCA_projPC{i,3}(hh,1,1),...
                'o', 'color',heading_colors(stimulus_type,:), 'MarkerSize',20, 'MarkerFaceColor',heading_colors(stimulus_type,:));
            % Start point of Choosing NULL
            h_null(hh) = plot3(PCA_projPC{i,1}(hh,2,1), PCA_projPC{i,2}(hh,2,1),PCA_projPC{i,3}(hh,2,1),...
                'o', 'color',heading_colors(stimulus_type,:), 'MarkerSize',20, 'linew',3);
            
            % PREF-choosing trajectory
            plot3(squeeze(PCA_projPC{i,1}(hh,1,:)), squeeze(PCA_projPC{i,2}(hh,1,:)),squeeze(PCA_projPC{i,3}(hh,1,:)),...
                '-', 'color',heading_colors(hh,:), 'linew',1.5);
            % Null-choosing trajectory
            plot3(squeeze(PCA_projPC{i,1}(hh,2,:)), squeeze(PCA_projPC{i,2}(hh,2,:)),squeeze(PCA_projPC{i,3}(hh,2,:)),...
                ':', 'color',heading_colors(hh,:), 'linew',1.5);
            
            for pp = 2:length(plotInd)
                plot3(PCA_projPC{i,1}(hh,1,plotInd(pp)), PCA_projPC{i,2}(hh,1,plotInd(pp)),PCA_projPC{i,3}(hh,1,plotInd(pp)),...
                    'o', 'color',heading_colors(hh,:), 'MarkerFaceColor',heading_colors(hh,:), 'linew',0.1, 'MarkerSize',15);
                
                plot3(PCA_projPC{i,1}(hh,2,plotInd(pp)), PCA_projPC{i,2}(hh,2,plotInd(pp)),PCA_projPC{i,3}(hh,2,plotInd(pp)),...
                    'o', 'color',heading_colors(hh,:), 'MarkerFaceColor','none', 'linew',2, 'MarkerSize',15);
            end
            
            grid on;
%             axis off;
            
        end
        end
    end
end


%% Optional: estimating "signal variance"

explVar_dpca_regu_signal = dpca_explainedVariance(firingRatesAverage, W_dpca_regu, V_dpca_regu, ...
    'combinedParams', combinedParams, ...
    'Cnoise', Cnoise, 'numOfTrials', trialNum);

% % Note how the pie chart changes relative to the previous figure.
% % That is because it is displaying percentages of (estimated) signal PSTH
% % variances, not total PSTH variances. See paper for more details.
%
% dpca_plot(firingRatesAverage, W_dpca_regu, V_dpca_regu, @dpca_plot_default, ...
%     'explainedVar', explVar_dpca_regu_signal, ...
%     'marginalizationNames', margNames, ...
%     'marginalizationColours', margColours, ...
%     'whichMarg', whichMarg_dpca_regu,                 ...
%     'time', time{1},                        ...
%     'timeEvents', timeEvents,               ...
%     'timeMarginalization', 3,           ...
%     'legendSubplot', 16);
% SetFigure(12);

%% Optional: decoding
if sort_id == ALL_CorrectCHOICE || sort_id ==All_Choice
    decodingClasses = {[(1:ST)' (1:ST)'], repmat([1:2], [ST 1]), [], [(1:ST)' (ST+(1:ST))']};
elseif sort_id == CORRECT_ANGLE  
    if ~isempty(stimulus_type)
        decodingClasses = {[(1:H)' (1:H)'], repmat([1:2], [H 1]), [], [(1:H)' (H+(1:H))']};
    else
        decodingClasses = {repmat((1:H)', [ST*2 1]), repmat(repmat([1:2], [H 1]), [1,ST]), repmat((1:ST), [H*2,1]), []};
    end
end

accuracy = dpca_classificationAccuracy(firingRatesAverage, firingRates, trialNum, ...
    'lambda', optimalLambda, ...
    'combinedParams', combinedParams, ...
    'decodingClasses', decodingClasses, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 100, ...        % increase to 100
    'filename', 'tmp_classification_accuracy.mat');

% dpca_classificationPlot(accuracy, [], [], [], decodingClasses)
% SetFigure(12);

accuracyShuffle = dpca_classificationShuffled(firingRates, trialNum, ...
    'lambda', optimalLambda, ...
    'combinedParams', combinedParams, ...
    'decodingClasses', decodingClasses, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 100, ...        % increase to 100
    'numShuffles', 100, ...  % increase to 100 (takes a lot of time)
    'filename', 'tmp_classification_accuracy.mat');

dpca_classificationPlot(accuracy, [], accuracyShuffle, [], decodingClasses)
SetFigure(12);

componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg_dpca_regu);

dpca_plot(firingRatesAverage, W_dpca_regu, V_dpca_regu, @dpca_plot_default, ...
    'explainedVar', explVar_dpca_regu_signal, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg_dpca_regu,                 ...
    'time', time{1},                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16,                ...
    'componentsSignif', componentsSignif);
SetFigure(12);


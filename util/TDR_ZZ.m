% Dimentional reduction using Targeted Dimensionality Reduction (TDR)
% PCA + Linear regression + QR decomposition
% Refer to Mante 2013 Nature

function TDR_ZZ(group_data, period, not_enough_reps_sess)

if nargin < 2
    period = 1;
    not_enough_reps_sess = [];
elseif nargin < 3
    not_enough_reps_sess = [];
end

% Pre-definition
%%%%%%% Order corresponding to "Sort_Id" in TEMPO_GUI processing %%%%
ALL_CorrectCHOICE = 1; CORRECT_ANGLE = 2; CHOICE_DIFFICULT = 3; OUTCOME = 4; WRONG_ANGLE = 5; CORRECTNWRONG_ANGLE = 6; All_Choice = 7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colors = [41 89 204; 248 28 83; 14 153 46]/255;
heading_colors = [0 0 0; 0.3137 0.1961 0.1249; 0.6274 0.3921 0.2497; 0.9363 0.5851 0.3726; 1.0000 0.7812 0.4975];
stimtype_name = {'Vestibular'; 'Visual'; 'Combined'};
% == Figure default
set(0, 'defaultFigureRenderer', 'painters')
set(0, 'defaultFigureRendererMode', 'Manual')

% Colors for each heading in each condition 
for k = 1:3
    colors_angles{k} = colormap(gray);
    colors_angles{k} = ones(size(colors_angles{k},1),3) - colors_angles{k} .* repmat([1 1 1]-colors(k,:),size(colors_angles{k},1),1);
    colors_angles{k} = colors_angles{k}(round(linspace(20,length(colors_angles{k}),5)),:);
end
set(0, 'defaultAxesColorOrder', cell2mat(colors_angles'))

group_data(not_enough_reps_sess)=[];
PSTH_temp = [group_data.PSTH];

PSTH = squeeze(PSTH_temp(period, CORRECT_ANGLE:7:end, :));   % only correct trial and stimulus period 

%% Concatenate the neural activities into a matrix (N_unit * (N_condition * T))
time_window ={ [0 1500], [-300 200], [0 800] };  % stimulus, around saccade, after feedback
ts = PSTH{1}.ts;

in_window = ts >= time_window{period}(1) & ts <= time_window{period}(2);  

%% 1. Linear regression
% Definition of variable 
% unique_heading = [0 -0 1 -1 2 -2 4 -4 8 -8];   % positve indicates choosing PREF direction
% head_norm = zscore(unique_heading); 
choice = repmat([1 -1], 1,5);    
modality = [1, -1];  % Positive indicates vestibular 

for n = 1:length(PSTH)
    unique_heading = unique(group_data(n).trialInfo(:,2));
    if length(unique_heading) == 8  % no 0 degree
        temp_headings = unique_heading;
    else  
        temp_headings = [unique_heading(5) -unique_heading(5) unique_heading(6) unique_heading(4) unique_heading(7) unique_heading(3) unique_heading(8) unique_heading(2) unique_heading(9) unique_heading(1)];
    end
    head_norm = zscore(temp_headings); 
    
    % Firing rate (Smoothed)
    PSTH_raw = [cell2mat(PSTH{n,1}.raw); cell2mat(PSTH{n,2}.raw)]; 
    PSTH_window = PSTH_raw(:,in_window); 
    PSTH_norm = reshape(zscore(PSTH_window(:)), size(PSTH_window,1), size(PSTH_window,2)); 
    
    % Construct Design Matrix
    for c = 1:2   % only individual-modality condition
        cond_num = cellfun(@(x) size(x,1), PSTH{n,c}.raw);
        
        headings_temp = []; choices_temp = []; 
        for cn = 1: length(cond_num)
%             headings_temp = head_norm(cn) * ones(cond_num(cn), 1); 
            headings_temp = [headings_temp; head_norm(cn) * ones(cond_num(cn), 1)];
            
%             choice_temp = choice(cn) * ones(cond_num(cn), 1); 
            choices_temp = [choices_temp; choice(cn) * ones(cond_num(cn), 1)]; 
        end
        headings{c,1} = headings_temp;
        choices{c,1} = choices_temp; 
        modalities{c,1} = ones(sum(cond_num),1) * modality(c); 
    end
    all_headings = cell2mat(headings);
    all_choices = cell2mat(choices);
    all_modalities = cell2mat(modalities);
    all_constant = ones(size(PSTH_norm,1),1); 
    
    design_mat = [all_headings all_choices all_modalities all_constant]; 
    
    % FR = beta1*headings + beta2*choices + beta3*modalities +beta4*constant 
    beta(n,:,:) = design_mat \ PSTH_norm; 
end

%% 2. PCA
for c = 1:3
    for n = 1:length(PSTH)
        % N_neuron * (10 conditions * Time)   (0(PREF) 0(Null) 1(PREF) 1(Null) 2(PREF) 2(Null) ... )
        PSTH_conca_temp{c}(n,:) = reshape(PSTH{n,c}.ys(:,in_window)', [], 1);   
        
    end
end
PSTH_conca = cell2mat(PSTH_conca_temp);   % N_neuron * (10 conditions * 3 stim_type * Time)

beta(sum(isnan(PSTH_conca),2)>0, :, :) = []; 
PSTH_conca(sum(isnan(PSTH_conca),2)>0, :) =[]; 

% Z-score of firing 
norm_PSTH_conca = zscore(PSTH_conca, 0, 2); 

% PCA
% Space spanned by individual modality and then project combined trials
% into this space
[weights_PCA, score, latent, ~, PCA_explained] = pca(norm_PSTH_conca(:,1:size(norm_PSTH_conca,2)/3*2)');   

%% 3. Regression subpace
denoised_dim=15; 
% Denoised matrix from PCA
D = weights_PCA(:, 1:denoised_dim) * weights_PCA(:, 1:denoised_dim)'; 

% Project regression weights into subspace spanned by D
for tt = 1:size(beta, 3)
    beta_pca(:,:,tt) = D * squeeze(beta(:,:,tt));
    
    % Norms of each regression vector at each time point
    for v = 1:3   % variables
        norm_vector(v,tt) = norm(squeeze(beta_pca(:,v,tt)));
    end
end

% Variable vector with maximum norm is defined as de-noised regression
% vector
[~,ind] = max(norm_vector, [], 2); 
beta_max = [beta_pca(:,1,ind(1))  beta_pca(:,2,ind(2)) beta_pca(:,3,ind(3))]; 

%%  Orthogonalization
[Q, R] = qr(beta_max); 
weight_TDR = Q(:,1:3); 

%% Projection 
time_bin = sum(in_window); 
heading_num = size(PSTH{n,c}.raw,1); 
modality_num = 3; 
project_norm_psth = reshape((weight_TDR' * norm_PSTH_conca)',...
    time_bin, heading_num, modality_num, size(weight_TDR,2));
project_norm_psth = permute(project_norm_psth, [2 1 3 4]);   % 10 headings, 150 time point, 3 conditions, 3 axis 

%% Plotting
plotInd = round(linspace(1,time_bin, 8)); 

% 1. 3-d trajectory
figure; hold on; 
for c = 1:size(project_norm_psth,3)  % stim_type
    for h = 1:size(project_norm_psth, 1)/2  % headings
        % Start point of choosing PREF
%         plot3(project_norm_psth(2*h-1,1,c,1), project_norm_psth(2*h-1,1,c,2),project_norm_psth(2*h-1,1,c,3),...
%             'o', 'color',colors_angles{c}(h,:), 'MarkerSize',20, 'MarkerFaceColor',colors_angles{c}(h,:));
%         % Start point of choosing NULL
%         plot3(project_norm_psth(2*h,1,c,1), project_norm_psth(2*h,1,c,2),project_norm_psth(2*h,1,c,3),...
%             'o', 'color',colors_angles{c}(h,:), 'MarkerSize',20, 'MarkerFaceColor',colors_angles{c}(h,:));
        plot3(project_norm_psth(2*h-1,1,c,1), project_norm_psth(2*h-1,1,c,2),project_norm_psth(2*h-1,1,c,3),...
            'o', 'color','k', 'MarkerSize',20, 'MarkerFaceColor','k', 'linew',2);
        plot3(project_norm_psth(2*h,1,c,1), project_norm_psth(2*h,1,c,2),project_norm_psth(2*h,1,c,3),...
            'o', 'color','k', 'MarkerSize',20, 'MarkerFaceColor','none', 'linew',2);
        
        % PREF-choosing trajectory
        plot3(project_norm_psth(2*h-1,:,c,1), project_norm_psth(2*h-1,:,c,2),project_norm_psth(2*h-1,:,c,3),...
            '-', 'color',colors_angles{c}(h,:), 'linew',3)
        % NULL-choosing trajectory
        plot3(project_norm_psth(2*h,:,c,1), project_norm_psth(2*h,:,c,2),project_norm_psth(2*h,:,c,3),...
            '--', 'color',colors_angles{c}(h,:), 'linew',3)
        
        for pp = plotInd
            plot3(project_norm_psth(2*h-1,pp,c,1), project_norm_psth(2*h-1,pp,c,2),project_norm_psth(2*h-1,pp,c,3),...
                'o', 'color',colors_angles{c}(h,:), 'MarkerFaceColor',colors_angles{c}(h,:), 'linew',0.1, 'MarkerSize',20);
            
            plot3(project_norm_psth(2*h,pp,c,1), project_norm_psth(2*h,pp,c,2),project_norm_psth(2*h,pp,c,3),...
                'o', 'color',colors_angles{c}(h,:), 'MarkerFaceColor','none', 'linew',2, 'MarkerSize',20);
        end

    end
end
SetFigure();

% 2. 2-d trajectory
% == Heading and choice
figure; hold on; 
for c = 1:size(project_norm_psth,3)  % stim_type
    for h = 1:size(project_norm_psth, 1)/2  % headings
        
        % Start point of choosing PREF
%         plot(project_norm_psth(2*h-1,1,c,1), project_norm_psth(2*h-1,1,c,2),...
%             'o', 'color',colors_angles{c}(h,:), 'MarkerSize',20, 'MarkerFaceColor',colors_angles{c}(h,:));
%         % Start point of choosing NULL
%         plot(project_norm_psth(2*h,1,c,1), project_norm_psth(2*h,1,c,2),...
%             'o', 'color',colors_angles{c}(h,:), 'MarkerSize',20, 'MarkerFaceColor',colors_angles{c}(h,:));
        plot(project_norm_psth(2*h-1,1,c,1), project_norm_psth(2*h-1,1,c,2),...
            'o', 'color','k', 'MarkerSize',20, 'MarkerFaceColor','k', 'linew',2);
        plot(project_norm_psth(2*h,1,c,1), project_norm_psth(2*h,1,c,2),...
            'o', 'color','k', 'MarkerSize',20, 'MarkerFaceColor','none', 'linew',2);
        
        % PREF-choosing trajectory
        plot(project_norm_psth(2*h-1,:,c,1), project_norm_psth(2*h-1,:,c,2),...
            '-', 'color',colors_angles{c}(h,:), 'linew',1.5)
        % NULL-choosing trajectory
        plot(project_norm_psth(2*h,:,c,1), project_norm_psth(2*h,:,c,2),...
            ':', 'color',colors_angles{c}(h,:), 'linew',1.5)
        
        for pp = plotInd
            plot(project_norm_psth(2*h-1,pp,c,1), project_norm_psth(2*h-1,pp,c,2),...
                'o', 'color',colors_angles{c}(h,:), 'MarkerFaceColor',colors_angles{c}(h,:), 'linew',0.1, 'MarkerSize',15);
            
            plot(project_norm_psth(2*h,pp,c,1), project_norm_psth(2*h,pp,c,2),...
                'o', 'color',colors_angles{c}(h,:), 'MarkerFaceColor','none', 'linew',2, 'MarkerSize',15);
        end
    end
end


% == Choice and modality
figure; hold on; 
for c = 1:size(project_norm_psth,3)  % stim_type
%     for h = 1:size(project_norm_psth, 1)/2  % headings
    for h = size(project_norm_psth, 1)/2  % headings
        
%         % Start point of choosing PREF
%         plot(project_norm_psth(2*h-1,1,c,2), project_norm_psth(2*h-1,1,c,3),...
%             'o', 'color',colors_angles{c}(h,:), 'MarkerSize',15, 'MarkerFaceColor',colors_angles{c}(h,:));
%         % Start point of choosing NULL
%         plot(project_norm_psth(2*h,1,c,2), project_norm_psth(2*h,1,c,3),...
%             'o', 'color',colors_angles{c}(h,:), 'MarkerSize',15, 'MarkerFaceColor','none');
        % Start point of choosing PREF
        plot(project_norm_psth(2*h-1,1,c,2), project_norm_psth(2*h-1,1,c,3),...
            'o', 'color','k', 'MarkerSize',20, 'MarkerFaceColor','k', 'linew',2);
        % Start point of choosing NULL
        plot(project_norm_psth(2*h,1,c,2), project_norm_psth(2*h,1,c,3),...
            'o', 'color','k', 'MarkerSize',20, 'MarkerFaceColor','none', 'linew',2);

        % PREF-choosing trajectory
        plot(project_norm_psth(2*h-1,:,c,2), project_norm_psth(2*h-1,:,c,3),...
            '-', 'color',colors_angles{c}(h,:), 'linew',1.5)
        % NULL-choosing trajectory
        plot(project_norm_psth(2*h,:,c,2), project_norm_psth(2*h,:,c,3),...
            ':', 'color',colors_angles{c}(h,:), 'linew',1.5)
        
        % -- Time markers
colorsHsv = repmat(rgb2hsv(colors(c,:)),length(plotInd),1);
colorsHsv(:,2) = linspace(0.2,1,length(plotInd));

if c == 3
    colorsHsv(:,3) = linspace(.9,colorsHsv(1,3),length(plotInd));
end

colorsRGB = hsv2rgb(colorsHsv);
colorsRGB = flipud(colorsRGB); 

        for pp = 1: length(plotInd)
            plot(project_norm_psth(2*h-1,plotInd(pp),c,2), project_norm_psth(2*h-1,plotInd(pp),c,3),...
                'o', 'color',colorsRGB(pp,:), 'MarkerFaceColor',colorsRGB(pp,:), 'linew',0.1, 'MarkerSize',15);
            
            plot(project_norm_psth(2*h,plotInd(pp),c,2), project_norm_psth(2*h,plotInd(pp),c,3),...
                'o', 'color',colorsRGB(pp,:), 'MarkerFaceColor','none', 'linew',2, 'MarkerSize',15);
        end
    end
end

% 3. 1-d trajectory
for ax = 1:3  
    figure; hold on; 

    for c = 1:size(project_norm_psth,3)  % stim_type
        for h = 1:size(project_norm_psth, 1)/2  % headings
                        
            % PREF-choosing trajectory
            plot(project_norm_psth(2*h-1,:,c,ax),...
                '-', 'color',colors_angles{c}(h,:), 'linew',1.5)
            % NULL-choosing trajectory
            plot(project_norm_psth(2*h,:,c,ax),...
                ':', 'color',colors_angles{c}(h,:), 'linew',1.5)
            
        end
    end
end

end

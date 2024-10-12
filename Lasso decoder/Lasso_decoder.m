% Training Lasso decoders across different time windows in different conditions
% Comsuming much time to train

% Input
% x, m*t*n dimension firing rates data
%     m cells simulateneouly recorded
%     t time windows, each time windows will be trained independently
%     n trials for training
% y, response variable, a vector with n categorical elements (0,1 or 1,-1)
% bootN, number of bootstraps, default is 100
% ifNorm, 0 or 1, if activity of each neuron has been zscored
% ifConstant , 0 or 1, if constant term has been included
% 'CV', n-fold cross-validation when training, default is 10
% variagin

% OutPut
% B, betas of each regressor (FR of each neuron)
% FitInfo, fitting information of each fitting, refer to function lasso document
function [B,FitInfo] = Lasso_decoder(x, y,bootN,ifNorm,ifConstant,CV,varagin)

nvaragin


if size(x,3) ~= length(y)
    error('Number of observation (trial number) is inconsistent');
end



find_for_Lasso = find(select_for_Lasso);
        
        for j = 1:length(j_for_Lasso)
            Lasso_tcenters{j} = min(rate_ts{j_for_Lasso(j)}) + Lasso_sliding_window/2 : Lasso_step_size : max(rate_ts{j_for_Lasso(j)}) - Lasso_sliding_window/2;
            
            % Regroup pseudo_trial_pool_decoder
            pseudo_trial_pool_decoder{j} = cell(length(select_for_Lasso),6,length(Lasso_tcenters{j}));
            
            for ii = 1:length(find_for_Lasso)
                
                raw_this = group_result(find_for_Lasso(ii)).mat_raw_PSTH.PSTH{j_for_Lasso(j),All_Choice,1}.raw;   % all trials
                
                for ttt = 1:length(Lasso_tcenters{j})
                    % Decoder windows
                    Lasso_epoch_ind = Lasso_tcenters{j}(ttt)-Lasso_sliding_window/2 <= rate_ts{j_for_Lasso(j)} ...
                        & rate_ts{j_for_Lasso(j)} <= Lasso_tcenters{j}(ttt)+Lasso_sliding_window/2;
                    
                    means_this = cellfun(@(x) mean(x(:,Lasso_epoch_ind),2), raw_this, 'UniformOutput', false);
                    
                    [pseudo_trial_pool_decoder{j}{ii, :, ttt}] = means_this{:};
                end
            end
            
            % Find cells who cross the minimal repetitions for each condition
            lengths = cellfun(@(x) length(x), pseudo_trial_pool_decoder{j}(:,:,1));
            select_for_decoder_decoder = all(lengths >= min_reps_decoder, 2);
            pseudo_trial_pool_decoder{j} = pseudo_trial_pool_decoder{j}(select_for_decoder_decoder,:, :);
            
            lengths_decoder = lengths(select_for_decoder_decoder,:);
            
            %----------------  Teacher signals (Response) ---------------------
            % Choice, 1 = Contralateral, 0 = Ipsilateral
            answers_choice = [ones(1,min_reps_training) zeros(1,min_reps_training)...
                ones(1,min_reps_training) zeros(1,min_reps_training) ones(1,min_reps_training) zeros(1,min_reps_training)]';

            % Get minimal reps in each conditions for all cells
            reps_actual_decoder = min(lengths_decoder);  
            
        end
        fprintf('Training Lasso decoder supposing that the %d neurons recorded simultaneously, \n', sum(select_for_decoder_decoder));
        fprintf('within which there are %d typical cells. \n' ,sum(select_for_decoder_decoder&select_tcells));
        
        % Initiation
        pseudo_trial_pool_perm = cell(length(j_for_Lasso),Lasso_bootstrap);
        training_set = cell(length(j_for_Lasso),Lasso_bootstrap);
        B = cell(length(Lasso_tcenters{j}),2,Lasso_bootstrap,size(pseudo_trial_pool_decoder{j},2)/2);
        FitInfo = cell(length(Lasso_tcenters{j}),2,Lasso_bootstrap,size(pseudo_trial_pool_decoder{j},2)/2);

        
        % Parallel computing
        if ver_num < 2014
            if matlabpool('size') == 0 ;  matlabpool;  end
        else
            if isempty(gcp('nocreate')); parpool; end
        end
        
        for j = 1:length(j_for_Lasso)
            progressbar('Lasso Bootstrap');
            for nn = 1:Lasso_bootstrap
                pseudo_trial_pool_perm{j,nn} = cellfun(@(x) x(randperm(size(x,1)),1), pseudo_trial_pool_decoder{j}, 'UniformOutput',false); % random permute trials in every condition and every time windows
                
                for conds = 1:size(pseudo_trial_pool_decoder{j},2)
                    training_set{j,nn}{conds} = reshape(cell2mat(squeeze(cellfun(@(x) x(1:min_reps_training,:,:), pseudo_trial_pool_perm{j,nn}(:,conds,:),'UniformOutput',false))),min_reps_training,[],length(Lasso_tcenters{j}));
                end
                
                for c = 1:size(pseudo_trial_pool_decoder{j},2)/2
                    firing_for_this_training = [training_set{j,nn}{2*c-1}; training_set{j,nn}{2*c}];
                    
                    parfor_progress(length(Lasso_tcenters{j}));
                    parfor ttt = 1:length(Lasso_tcenters{j})
                        % z-score for every cell
                        norm_trainFR = zscore(firing_for_this_training(:,:,ttt),0,2);
                        
                        firing_for_training = [norm_trainFR ones(size(firing_for_this_training,1),1)];  % Add a constant term
                        
                        % LASSO Fitting
                        %                     options = statset('UseParallel',1);
                        rng default;  % For reproducibility
                        [B{ttt,j,nn,c},FitInfo{ttt,j,nn,c}] = lasso(firing_for_training,answers_choice((c-1)*80+1:c*80), 'CV', 10);
                        
                        parfor_progress;
                    end
                    parfor_progress(0);
                end
                progressbar(nn/Lasso_bootstrap);
            end
        end

%% ================= Testing  Part  =====================
% norm_testFR = cell(length(Lasso_tcenters),Lasso_bootstrap,6);

for j = 1:length(Lasso_tcenters)
    for nn = 1: Lasso_bootstrap
        
        for c = 1:6
            % Inexperienced trials
            testing_set{j,nn,c}= reshape(cell2mat(squeeze(cellfun(@(x) x(min_reps_training+1:reps_actual_decoder(c),:,:),...
                pseudo_trial_pool_perm{j,nn}(:,c,:),'UniformOutput',false))),...                            % cell2mat
                (reps_actual_decoder(c)-min_reps_training),[],length(Lasso_tcenters{j}));       % reshape
            
            % z-score
            norm_testFR{j,nn,c} = zscore(testing_set{j,nn,c},0,2);
        end
        
        for pc = 1:3   % paired condition
            % Each two conditions share a same decoder (same decoder in the same stim type)
            
            fitinfo_temp = FitInfo(1:length(Lasso_tcenters{j}),j,nn,pc);
            minMSE_ind = squeeze(cellfun(@(x) x.IndexMinMSE, fitinfo_temp));   % Min MSE of 10-fold cross-validation test
        
            for ttt = 1:length(Lasso_tcenters{j})
                % Capture the Betas with minimal MSE in 10-fold CV training
                B_MinMSE{ttt,j,nn,pc} = B{ttt,j,nn,pc}(:,minMSE_ind(ttt));   
                
                % Add the constant terms
                norm_testFR_constInclud{j,nn,pc}(:,:,ttt) = [norm_testFR{j,nn,2*pc-1}(:,:,ttt) ones(size(testing_set{j,nn,2*pc-1},1),1); ...
                    norm_testFR{j,nn,2*pc}(:,:,ttt) ones(size(testing_set{j,nn,2*pc},1),1)];
                
                % Decode the choice window-by-window from the normalized FR
                predict_choice{j,nn,pc}(:,ttt) = norm_testFR_constInclud{j,nn,pc}(:,:,ttt) * B_MinMSE{ttt,j,nn,pc};
            end
            
            predict_contral{j,nn,pc} = predict_choice{j,nn,pc}(1:size(testing_set{j,nn,2*pc-1},1),:) >= 0.5;  % Correctly predict contrallater choice
            predict_ipsi{j,nn,pc} = predict_choice{j,nn,pc}(size(testing_set{j,nn,2*pc-1},1)+1:end,:) <= 0.5; % Correctly predict ipsilateral choice
            
            % Prediction correct rate in each conditions
            predict_accuracy{j,nn,pc} = (sum(predict_contral{j,nn,pc})+sum(predict_ipsi{j,nn,pc}))/size(predict_choice{j,nn,pc},1);
            
        end
        
    end
end
        
    

for j =1:length(Lasso_tcenters)
    for pc = 1:size(pseudo_trial_pool_decoder{1},2)/2
       grouped_pred_accuracy{j,pc} = reshape(cell2mat(squeeze(predict_accuracy(j,:,pc))),[], nn)';
       mean_pred_accu{j,pc} = mean(grouped_pred_accuracy{j,pc});
       std_pred_accu{j,pc} = std(grouped_pred_accuracy{j,pc});
       
       group_B_MinMSE{j,pc} = reshape(cell2mat(squeeze(B_MinMSE(1:length(Lasso_tcenters{j}),j,:,pc))),[],length(Lasso_tcenters{j}),nn);
       mean_B_MinMSE{j,pc} = mean(group_B_MinMSE{j,pc},3);
       std_B_MinMSE{j,pc} = std(group_B_MinMSE{j,pc},0,3);
       
    end
end

%% Saving
% Because the huge cost of time I saving the variables for future
cd
if exist('Lasso_Decoder.mat')
    delete 'Lasso_Decoder.mat' ;
end
save('Lasso_Decoder.mat', 'grouped_pred_accuracy','mean_pred_accu','std_pred_accu',...
    'group_B_MinMSE','mean_B_MinMSE','std_B_MinMSE',...
    'B','FitInfo');


%% ========== Plotting ================
% --- 1. Predicted accuracy
% Visualizing prediction accuracy of each bootstrap
set(figure(figN), 'pos',[10 20 1200, 900], 'Name', 'Prediction Accuracy of Each Reps.'); clf; figN=figN+1;
ratio_w = [length(Lasso_tcenters{1})  length(Lasso_tcenters{2})];
ha = tight_subplot(3,2,[.05 .01], [.1 .1], [.1 .05],[],ratio_w);
for pc = 1:size(pseudo_trial_pool_decoder{1},2)/2
    for j = 1:length(Lasso_tcenters)
        set(gcf,'CurrentAxes', ha(pc+3*(j-1)));hold on;
        
%         clims = [0 1];
        imagesc(grouped_pred_accuracy{j,pc},[0 1]);
        
        if pc==1 && j==2
            colorbar('northoutside');
        end
    end
end
 





% Averaged across 1000 bootstraps
set(figure(figN), 'pos',[10 20 1600, 900], 'Name', 'Prediction Accuracy of Lasso Decoder'); clf; figN=figN+1;
ratio_w = [length(Lasso_tcenters{1})  length(Lasso_tcenters{2})];
ha = tight_subplot(3,2,[.05 .01], [.1 .1], [.1 .05],[],ratio_w);
for pc = 1:size(pseudo_trial_pool_decoder{1},2)/2
    for j = 1:length(Lasso_tcenters)
        set(gcf,'CurrentAxes', ha(pc+3*(j-1)));hold on;
        plot(mean_pred_accu{j,pc},'Color',colors(pc,:),'LineWidth',2);
        plot(mean_pred_accu{j,pc}+std_pred_accu{j,pc},'--', 'Color',colors(pc,:),'LineWidth',0.5);
        plot(mean_pred_accu{j,pc}-std_pred_accu{j,pc},'--', 'Color',colors(pc,:),'LineWidth',0.5);

    end
end

end
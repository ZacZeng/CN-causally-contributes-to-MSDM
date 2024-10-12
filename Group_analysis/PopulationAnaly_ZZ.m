% Neural Population Analysis
% Aim to population data from any Batch directory

% 1. Mean FR
% 2. dPCA
% 3. LASSO decoder
% 4. SVM

function PopulationAnaly_ZZ(~)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======= Change the file directory===========

cd 'D:\Paper_rawdata\Raw_data\LIP\LIP_m5_m10';
% cd 'D:\Paper_rawdata\Raw_data\FEF\FEF_m10_m13';
% cd 'D:\Paper_rawdata\Raw_data\MSTd\MSTd_GY_2008';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Order corresponding to "Sort_Id" in TEMPO_GUI processing %%%%
ALL_CorrectCHOICE = 1; CORRECT_ANGLE = 2; CHOICE_DIFFICULT = 3; OUTCOME = 4; WRONG_ANGLE = 5; CORRECTNWRONG_ANGLE = 6; All_Choice = 7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Pool all sessions
fileName_list = dir('*_PSTH.mat');
not_enough_reps_sess = [];

progressbar('Load .mat files');
for sess = 1:size(fileName_list,1)
    load(fileName_list(sess).name);
    
    if result.repetitionN < 8 || length(result.unique_stim_type) ~= 3 || length(unique(result.heading_per_trial)) < 9 % Change to 8 when run FEF data £¡£¡£¡£¡£¡£¡£¡£¡£¡£¡£¡
        not_enough_reps_sess = [not_enough_reps_sess ; sess];
    end
    
    group_result(sess).file = [result.FILE '_' num2str(result.SpikeChan)];
    group_result(sess).stimtype = result.unique_stim_type;
    group_result(sess).trialInfo = [result.stim_type_per_trial' result.heading_per_trial' result.choice_per_trial result.outcome_per_trial];
    group_result(sess).pref = result.PREF;
    group_result(sess).rawSpk = result.spike_aligned;
    group_result(sess).align_markers = result.align_markers;
    group_result(sess).align_offsets_others = result.align_offsets_others;
    group_result(sess).PSTH = result.PSTH;
    
    group_result(sess).choiceDive = result.ChoiceDivergence_ALL;
    group_result(sess).choiceDive_p = result.ChoiceDivergence_ALL_perm;
    group_result(sess).choiceDive_Diffic = result.ChoiceDivergence_Difficult;
    group_result(sess).choiceDive_Easy = result.ChoiceDivergence_Easy;
    group_result(sess).choiceDive_Correct = result.ChoiceDivergence_Correct;
    group_result(sess).choiceDive_Error = result.ChoiceDivergence_Error;
    group_result(sess).choiceDive_CorrectHeading = result.ChoiceDivergence_CorrectHeading;
    group_result(sess).choiceDive_ErrorHeading = result.ChoiceDivergence_ErrorHeading;
    
    group_result(sess).choicePref = result.ChoicePreference;
    group_result(sess).choicePref_p = result.ChoicePreference_pvalue;
    group_result(sess).choicePref_correctangle = result.ChoicePref_correctangle;
    group_result(sess).choicePref_correctangle_p = result.ChoicePref_correctangle_p;
    group_result(sess).choicePref_errorangle = result.ChoicePref_errorangle;
    group_result(sess).choicePref_errorangle_p = result.ChoicePref_errorangle_p;   % Error spelling hhh
    group_result(sess).choicePref_difficeasy = result.ChoicePref_difficeasy;
    group_result(sess).choicePref_difficeasy_p = result.ChoicePref_difficeasy_p;
    group_result(sess).choicePref_correcterror = result.ChoicePref_correcterror;
    group_result(sess).choicePref_correcterror_p = result.ChoicePref_correcterror_p;
    
    group_result(sess).modalityPref = result.ModalityPreference; 
    
    group_result(sess).CP = result.CP;
    
    % For Messi Reaction-time task
    % Added by ZZ @20230518
    m_ind = find(result.FILE == 'm'); 
    c_ind = find(result.FILE == 'c'); 
    r_ind = find(result.FILE == 'r'); 
    if str2num(result.FILE(m_ind+1 : c_ind-1)) == 10 && str2num(result.FILE(c_ind+1 : r_ind-1)) > 1500
        group_result(sess).RT_median_condition = result.RT_median_condition; 
        group_result(sess).RT_median_heading = result.RT_median_heading; 
        group_result(sess).unique_heading = result.RT_unique_heading; 
    end
    
    progressbar(sess/size(fileName_list,1));
end


%% Dimensional reduction
% PCA, dPCA
% for st = 1:3
%     %         dPCA_ZZ(group_result, CORRECTNWRONG_ANGLE,st);
%     dPCA_ZZ(group_result, CORRECT_ANGLE,st);
% end
dPCA_ZZ(group_result, ALL_CorrectCHOICE,[]);

% % Gouki's suggestion: trajectory of three modalities in the same dPCA subspace
% dPCA_ZZ(group_result, CORRECT_ANGLE,[]);

% TDR
period = 1;  % 1, stimulus period; 2, around sccade; 3, after feedback
TDR_ZZ(group_result, period, not_enough_reps_sess); 

%% Decoding
% SVM or Lasso
% ---------------------------------------------- Data preparing for decoder -------------------------------------------------
if exist('decoder_prepare_data.mat')
    load('decoder_prepare_data.mat');
else
    % Data preparing for Decoder training
    decoder_data = decoder_training_prepare(group_result,All_Choice,'bootstrapN', 100);  % Changed from 1000 to 100 for the sake of memory capacity
    
    clear group_result;   % Free memory
    
%     save('decoder_prepare_data.mat', 'decoder_data', '-v7.3');
end


% ---------------------------------------------- Training -------------------------------------------------------

if exist('trained_decoder.mat')
    load('trained_decoder.mat');
else
    % Initiation
    model= cell(2,3);
    SVM = 1;   LASSO_A = 2;
    
%     for type = SVM:LASSO_A
    for type = LASSO_A
        for st = 1:3
            model{type, st} = decoder_train(decoder_data.FR,decoder_data.training_set, decoder_data.teacher_signal,st,type);
        end
    end
    
%     save('trained_decoder.mat', 'model', '-v7.3');
end


%  Plot correlation of each decoder in different windows
SVM = 1;   LASSO_A = 2;
% load('decoder_tcenters.mat');

progressbar('Plot Decoder Weight')
% for type = SVM : LASSO_A
for type = LASSO_A
    for st = 1:3
        plot_decoder_correlation(model{type,st}, type, st, decoder_data.t_centers);
        
        progressbar(((type-1)*3+st) / 6);
    end
end


%% Across-modality decoding
% only for Lasso
% load('decoder_prepare_data.mat');
% load('trained_decoder.mat');
% load('decoder_tcenters.mat');
type = 2; 
for  st_train = 1:3
    for st_test = 1:3
        accuracy_training{st_train, st_test} = decoder_test(model{type,st_train},decoder_data.FR, decoder_data.training_set,st_test,type);
        
        progressbar(((st_train-1)*3 +st_test)/9);
    end
end

for st_train = 1:3
    for st_test = 1:3
                plot_temporal_decoder_accuracy(accuracy_training{st_train,st_test}, type, st_test, decoder_data.t_centers);
    end
end


% --------------------------------------------------------- Testing -----------------------------------------------
% Predicting accuracy for testing set
% Cross temproal decoding using training set
% if exist('accuracy_test_set.mat') && exist('accuracy_training_set.mat')
%     load('accuracy_test_set.mat');
%     load('accuracy_training_set.mat');
%     
% else
    SVM = 1;   LASSO_A = 2;
    
    for type = LASSO_A
        progressbar('Decoder Testing')
        for st = 1:3
%             accuracy_training{type, st} = decoder_test(model{type,st},decoder_data.training_set,st,type);
%             accuracy_test{type, st} = decoder_test(model{type,st},decoder_data.testing_set,st,type);
                accuracy_test{type, st} = decoder_test(model{type,st},decoder_data.FR, decoder_data.testing_set,st,type);
            progressbar(st/3);
        end
    end
    
%     save('accuracy_test_set.mat', 'accuracy_test', '-v7.3');
%     save('accuracy_training_set.mat', 'accuracy', '-v7.3');
% end


% Plot cross-temporal decoder accuracy
% load('decoder_tcenters.mat');
SVM = 1;   LASSO_A = 2;

progressbar('Plot Performance of Decoder')
for type = LASSO_A
    for st = 1:3
%         plot_temporal_decoder_accuracy(accuracy_training{type,st}, type, st, decoder_tcenters);
                plot_temporal_decoder_accuracy(accuracy_test{type,st}, type, st,decoder_data.t_centers);
        
        progressbar(((type-1)*3+st) / 6);
    end
end


end

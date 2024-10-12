% Linear decoder training (SVM or Lasso)

% FR, firing rates for decoders
% Added by ZZ @ 20231121

% xx, training set of predictor, output from script of
% decoder_training_prepare.m. J*nn cell, j is the way of alignment and nn
% is the number of bootstrap
% changed from FR to only trial index @ 20231121

% y, answers /  response. Categorical values.

% st, stimulus type. Just as x, y includes all 3 stimulus data. So here, I
% indicate which stimulus type I want to train, so that I can extract data
% from corresponding types.

% decoder_type, SVM (1) or Lasso (2) now.
function model = decoder_train(FR,xx,y,st,decoder_type)

% If you come from script decoder_training_prepare.m, and cast the
% decoder_data.training_set into this script, then we can now the j,
% bootstrapN, trial repetitions from the x.
j_for_decoder = size(xx,1);
bootstrapN = size(xx,2);
reps = size(xx{j_for_decoder,bootstrapN}{1,1},2);
% reps = size(x{j_for_decoder,bootstrapN}{1,st,1},1);


if decoder_type == 1   % SVM
    % Initiation
    SVMModel = cell( j_for_decoder,bootstrapN,size(xx{1,1},3));
    % Pre-defining for optimization setting
    opts = struct('Optimizer', 'bayesopt','Kfold',10,'MaxObjectiveEvaluations',30,...
        'ShowPlots',0, 'AcquisitionFunctionName','expected-improvement-plus');
    
    progressbar('Trainning Bootstrap');
    for nn = 1:bootstrapN
        for j = 1:j_for_decoder
            predictor = xx{j,nn};
            
            parfor_progress(size(predictor,3));
            parfor tt = 1: size(predictor,3)
                
                %                 rng default;  % For reproducibility
                %                 SVMModel{j,nn,tt} = fitcsvm(cell2mat(predictor(:,2*st-1:2*st,tt)'), y(2*reps*(st-1)+1:2*reps*st),...
                %                     'Standardize',true,'OptimizeHyperparameters','auto',...
                %                     'HyperparameterOptimizationOptions', opts);
                
                % Delete cross-validation and regularization for extremely killing time
                SVMModel{j,nn,tt} = fitcsvm(cell2mat(predictor(:,2*st-1:2*st,tt)'), y(2*reps*(st-1)+1:2*reps*st),...
                    'Standardize',true);
                
                
                parfor_progress;
            end
            parfor_progress(0);
        end
        progressbar(nn/bootstrapN);
        
    end
    
    model = SVMModel;
    
elseif decoder_type == 2  % Lasso
    % Referring to Kiani et al., 2014, CB
    
    % Initiation
    B = cell( j_for_decoder,bootstrapN,size(FR{1,1},3));
    %     B = cell( j_for_decoder,bootstrapN,size(x{1,1},3));
    FitInfo = B;
    
    %     progressbar('Lasso Bootstrap');
    h = waitbar(0, ' ', 'Name', 'Lasso is training');
   
    
    %     for nn = 1:bootstrapN
    for j =1: j_for_decoder
        
        %             predictor = x{j,nn};
        
        %             parfor_progress(size(predictor,3));
        %             parfor tt = 1:size(predictor,3)
        for tt = 1:size(FR{j},3)
            
                % Where is the decoder
                msg = [' j = ', num2str(j), ' tt = ', num2str(tt)];
                if j==1
                    locat = j * tt;
                elseif j==2
                    locat = (j-1) * tt + size(FR{1},3);
                end
                waitbar(locat/ (size(FR{1},3)+size(FR{2},3)), h, msg);
            
                parfor_progress(bootstrapN);
            parfor nn = 1:bootstrapN
                
                predictor = cellfun(@(x, y) x(y), FR{j}(:,:,tt), xx{j,nn}, 'UniformOutput', 0);
                
                rng default;  % For reproducibility
                [bb, FitInfo] = lasso(cell2mat(predictor(:,2*st-1:2*st)'),y(2*reps*(st-1)+1:2*reps*st), 'CV', 10, 'Standardize',1);
                %                 [B{j,nn,tt}, FitInfo{j,nn,tt}] = lasso(cell2mat(predictor(:,2*st-1:2*st,tt)'),y(2*reps*(st-1)+1:2*reps*st), 'CV', 10, 'Standardize',1);
                
                min_mse_ind = FitInfo.IndexMinMSE;
                B{j,nn,tt} = bb(:,min_mse_ind);
                intercep{j,nn,tt} = FitInfo.Intercept(min_mse_ind);
                                
                                parfor_progress;
            end
                        parfor_progress(0);
            
        end
        %         progressbar(nn/bootstrapN);
        
    end
    
    delete(h);

    model.B= B;
    model.intercep = intercep;
    %     for j = 1: j_for_decoder
    %         for nn = 1: bootstrapN
    %             predictor = xx{j,nn};
    %             for tt = 1: size(predictor,3)
    %
    %                 min_mse_ind = FitInfo{j,nn,tt}.IndexMinMSE;
    %                 model.B{j,nn,tt} = B{j,nn,tt}(:,min_mse_ind);
    %                 model.intercep{j,nn,tt} = FitInfo{j,nn,tt}.Intercept(min_mse_ind);
    %
    %             end
    %         end
    %     end
    
    %     model.B = B;
    %     model.FitInfo = FitInfo;
end
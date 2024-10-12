% Linear decoder testing (SVM or Lasso)

% model, Models of SVM or Lasso.

% FR, FR of each temporal window

% xx, testing set of predictor, output from script of
% decoder_training_prepare.m. J*nn cell, j is the way of alignment and nn
% is the number of bootstrap

% % % % % % y, answers /  response. Categorical values. 6 cells with cell (2*st-1) of them are
% % % % % % positive trials while the Num. 2*st are negative ones.

% st, stimulus type. Just as x, y includes all 3 stimulus data. So here, I
% indicate which stimulus type I want to train, so that I can extract data
% from corresponding types.

% decoder_type, SVM (1) or Lasso (2) now.
function correct_rate = decoder_test(model,FR, xx,st,decoder_type)

% If you come from script decoder_training_prepare.m, and cast the
% decoder_data.training_set into this script, then we can now the j,
% bootstrapN, trial repetitions from the x (reps are different in different conditions).
j_for_decoder = size(xx,1);
bootstrapN = size(xx,2);
% reps = size(x{j_for_decoder,bootstrapN}{1,st,1},1);


if decoder_type == 1  % SVM
    
    h = waitbar(0, ' ', 'Name', 'SVM is decoding');
    
    %     fig = uifigure;
    %     d = uiprogressdlg(fig, 'Title', ' SVM is decoding ', 'Message', '1', 'Cancelable','on');
    %     drawnow
    
    %     parfor_progress(bootstrapN);
    for nn = 1:bootstrapN
        
        %         % progessbar massage
        %         d.Value = nn / bootstrapN;
        
        for j = 1:j_for_decoder
            for tt = 1: size(xx{j,nn},3)
                model_svm = model{j,nn,tt};
                
                for ttt = 1:size(xx{j,nn},3)    % Cross-temporal decoding
                    
                    %                     % where the decoder is
                    msg = ['bootN = ' num2str(nn) ' j = ' num2str(j) '  trainT = ' num2str(tt) ' testT = ' num2str(ttt)];
                    waitbar(nn / bootstrapN, h, msg);
                    %                     msg = ['j = ' num2str(j) '  trainT = ' num2str(tt) ' testT = ' num2str(ttt)];
                    %                     d.Message = msg;
                    
                    predictor = [xx{j,nn}{:,st*2-1,ttt} ; xx{j,nn}{:,st*2,ttt} ];     % testing set of other time windows
                    
                    oneN = size(xx{j,nn}{1,st*2-1,ttt},1);     % Number of trials should answer 1
                    zeroN = size(xx{j,nn}{1,st*2,ttt},1);        % Number of trials should answer 0
                    
                    %                     response = [y{st*2-1}(:) ; y{st*2}(:)];
                    response = [ones(oneN,1); zeros(zeroN,1)];
                    
                    correct_rate{j,nn}(tt,ttt) = 1 - loss(model_svm, predictor, response);
                end
            end
        end
        %         parfor_progress;
        
    end
    %     parfor_progress(0);
    %     close(d);
    delete(h);
    
elseif decoder_type == 2   % Lasso
    
    h = waitbar(0, ' ', 'Name', 'Lasso is decoding');
    %     fig = uifigure;
    %     d = uiprogressdlg(fig, 'Title', ' Lasso is decoding ', 'Message', '1','Cancelable','on');  % supported only after 2018a
    %     drawnow
    
    
    %         % progessbar massage
    %         d.Value = nn / bootstrapN;    
    for j = 1:j_for_decoder
        for tt = 1: size(FR{j},3)
            
            % where the decoder is
            msg = [' j = ' num2str(j) '  trainT = ' num2str(tt)];
            waitbar(j*tt / (j * size(FR{j},3)) , h, msg);
            %                     d.Message = msg;
            
            
            %             for tt = 1: size(xx{j,nn},3)
            % Extract minimal cross-validation model of the current time window
            %                 min_error_ind = model.FitInfo{j,nn,tt}.IndexMinMSE;               % minimal cross-validation error trial index
            %                 intercep =  model.FitInfo{j,nn,tt}.Intercept(min_error_ind);    % intercept
            %                 beta_minMSE = model.B{j,nn,tt}(:,min_error_ind);               % coefficients of each predictors in the best trial
                
%                 for nn = 1:bootstrapN
                    
                    % Change for loop to martix multiplication
                    % ZZ @20231129
                    temp_intercep = [model.intercep{j,:,tt};];
                    temp_beta = [model.B{j,:,tt}];
                    
                    xxx = xx(j,:); frs = FR{j}; 
                    correct_rate_temp = []; 
                    
            parfor_progress(size(FR{j},3));       
            parfor ttt = 1:size(frs,3)   % Cross-temporal decoding
                    % It's a hard way, hope it deserves
                    
                    % Preferred choice
                    xxx_temp = cellfun(@(x) [x(:, st*2-1);], xxx, 'UniformOutput', 0);  % trial index
                    xxx_temp = [xxx_temp{:}]; 
                    FR_temp = repmat(squeeze(frs(:, st*2-1, ttt)),1, bootstrapN);
                    FR_temp = cell2mat(cellfun(@(x,y) x(y), FR_temp, xxx_temp, 'UniformOutput', 0));  % FR of needed trials
                    FR_temp = permute(reshape(FR_temp, length(xxx_temp{1,1}), size(xxx_temp,1), []), [2 1 3]); 
                    beta = reshape(repmat(temp_beta, size(FR_temp,2), 1), size(FR_temp));
                    answer_pref = squeeze(dot(FR_temp, beta)) + temp_intercep;
                    correct_pref = sum (answer_pref >= 0.5) / size(answer_pref,1);
                    
                    % NULL-Preferred choice
                    xxx_temp = cellfun(@(x) [x(:, st*2);], xxx, 'UniformOutput', 0);  % trial index
                    xxx_temp = [xxx_temp{:}]; 
                    FR_temp= repmat(squeeze(frs(:, st*2, ttt)),1, bootstrapN);
                    FR_temp = cell2mat(cellfun(@(x,y) x(y), FR_temp, xxx_temp, 'UniformOutput', 0));  % FR of needed trials
                    FR_temp = permute(reshape(FR_temp, length(xxx_temp{1,1}), size(xxx_temp,1), []), [2 1 3]); 
                    beta = reshape(repmat(temp_beta, size(FR_temp,2), 1), size(FR_temp));
                    answer_null = squeeze(dot(FR_temp, beta)) + temp_intercep;
                    correct_null = sum(answer_null < 0.5) / size(answer_null, 1);
                    
                    correct_rate_temp(ttt,:) = (correct_pref+correct_null)/2; 
                    
                    parfor_progress;

%                 end
            end
            correct_rate{j}(tt,:,:) = correct_rate_temp;
%             for ttt = 1:size(FR{j},3)   % Cross-temporal decoding
%                 
%                     % where the decoder is
%                     msg = [' j = ' num2str(j) '  trainT = ' num2str(tt) ' testT = ' num2str(ttt)];
%                     waitbar(j*tt*ttt / (j * size(FR{j},3) * size(FR{j},3)) , h, msg);
%                     %                     d.Message = msg;
%                 
% %                     parfor_progress(bootstrapN);
%                 for nn = 1:bootstrapN
%                     
%                     % Changed by ZZ @ 20231120
%                     intercep = model.intercep{j,nn,tt};
%                     beta_minMSE = model.B{j,nn,tt};
% 
%                     %                     predictor = [xx{j,nn}{:,st*2-1,ttt} ; xx{j,nn}{:,st*2,ttt} ];      % testing set of other time windows
%                     % added by ZZ @ 20231127
%                     predictor = ([ [xx{j,nn}(:,st*2-1);]  [xx{j,nn}(:,st*2);] ]);      % testing set of other time windows
%                     predictor = cellfun(@(x, y) x(y), FR{j}(:,st*2-1:st*2,ttt), predictor, 'UniformOutput', 0);
%                     predictor = cell2mat(predictor');
%                     
%                     oneN = size(xx{j,nn}{1,st*2-1},2);     % Number of trials should answer 1
%                     zeroN = size(xx{j,nn}{1,st*2},2);        % Number of trials should answer 0
%                     %                     oneN = size(xx{j,nn}{1,st*2-1,ttt},1);     % Number of trials should answer 1
%                     %                     zeroN = size(xx{j,nn}{1,st*2,ttt},1);        % Number of trials should answer 0
%                     
%                     %         y_predict(t,tt,:) = [[decoder_data.training_set{1}{:,1,tt} ; [decoder_data.training_set{1}{:,2,tt}]] ones(80,1)] * beta_minMSE(:,t) + intercep(:,t);
%                     y_predict = predictor* beta_minMSE+ intercep;
%                     
%                     corrctN_pref = sum(y_predict(1:oneN) >= 0.5) ;     % Prediction correct trial number in choosing PREF trials
%                     corrctN_null = sum(y_predict(oneN+1:end) <0.5);  % Prediction correct trial number in chossing NULL trials
% %                     correct_rate{j,nn}(tt,ttt) = (corrctN_pref+corrctN_null) / length(y_predict);
%                     
%                 parfor_progress;
%                 end
%             end
        end
    end
        parfor_progress(0);
    %     close(d);
    delete(h);
end

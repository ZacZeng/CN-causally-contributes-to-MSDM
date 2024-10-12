% Data preparing for Decoders (e.g. SVM, Lasso)
% @20231121
% Change: using the same trials for different temporal decoders in the same bootstraps

function decoder_prepare = decoder_training_prepare(X, sort_id, varargin)

% X is the group_data extracted from a Batch directory. It include many
% structures, each of which is from a cell. But we only need only a few
% field of them. sort_id tells PSTH data sorted by what is required. 

% varargins
% --  j_for_decoder, how many period will be analysed,
% from 1 to 3, default is 2

% -- min_reps4training, how many trials are need for decoder training in each
% conditions

% -- min_reps4each_condition, how many trials a cell need conclude in all
% conditions. If a cell doesn't have enough trials, we give up it. 

% -- decoder_window, how long period response we use to train a decoder

% -- decoder_step_size, intervals between two decoders 

% -- bootstrapN, number of bootstraps, default is100; 

% default input parameters
options = struct( 'j_for_decoder',      2, ...
                 'min_reps4training',     40, ...
                 'min_reps4each_condition',          50, ...
                 'decoder_window', 100, ...
                 'decoder_step_size', 20,...
                 'bootstrapN', 100);

% read input parameters
optionNames = fieldnames(options);
if mod(length(varargin),2) == 1
	error('Please provide propertyName/propertyValue pairs')
end
for pair = reshape(varargin,2,[])    % pair is {propName; propValue}
	if any(strcmp(pair{1}, optionNames))
        options.(pair{1}) = pair{2};
    else
        error('%s is not a recognized parameter name', pair{1})
	end
end


%% Extract PSTH
PSTH_temp = [X(:).PSTH];
PSTH_sort_id =PSTH_temp(1: options.j_for_decoder,sort_id:7:end,1);

rate_ts = {X(1).PSTH{1,1,1}.ts,  X(1).PSTH{2,1,1}.ts, X(1).PSTH{3,1,1}.ts};

clear 'X'; 

for j = 1: options.j_for_decoder
    decoder_tcenters{j} = min(rate_ts{j}) + options.decoder_window/2 : options.decoder_step_size : max(rate_ts{j}) - options.decoder_window/2;
    pseudo_trial_decoder{j} = cell(length(PSTH_sort_id),6,length(decoder_tcenters{j}));
end


for j = 1:options.j_for_decoder
    for ii = 1:size(PSTH_sort_id,2)
        raw_this = PSTH_sort_id{j,ii}.raw;
%         raw_this = cellfun(@(x) zscore(x,0, 2), raw_this, 'UniformOutput', 0);  % added 20240321
        
        for tt = 1:length(decoder_tcenters{j})
            % Decoder windows
            decoder_epoch_ind = decoder_tcenters{j}(tt)-options.decoder_window/2 <= rate_ts{j} ...
                & rate_ts{j} <= decoder_tcenters{j}(tt)+options.decoder_window/2;
            
            mean_this_window = cellfun(@(x) mean(x(:,decoder_epoch_ind),2), raw_this, 'UniformOutput', 0);
            
            [pseudo_trial_decoder{j}{ii,:,tt}] = mean_this_window{:};
        end
    end
    
    % Find cells who cross the minimal repetitions threshold in all condition
    reps = cellfun(@length, pseudo_trial_decoder{j}(:,:,1));
    enough_reps_cell = all(reps >= options.min_reps4each_condition, 2);
    
    pseudo_trial_decoder{j} = pseudo_trial_decoder{j}(enough_reps_cell,:,:);
    reps_decoder = reps(enough_reps_cell,:);
    
end

%----------------  Teacher signals (Response) ---------------------
% Choice, 1 = Contralateral, 0 = Ipsilateral
answers_choice = [ones(1,options.min_reps4training) zeros(1,options.min_reps4training)...
    ones(1,options.min_reps4training) zeros(1,options.min_reps4training) ones(1,options.min_reps4training) zeros(1,options.min_reps4training)]';

% Get minimal reps in each conditions for all cells
min_reps_each_condition = min(reps_decoder);

fprintf('Training Lasso decoder supposing that the %d neurons recorded simultaneously, \n',sum(enough_reps_cell));


%% == Bootstraps ==
% Clear some variables to free memory
clear X PSTH_temp PSTH_sort_id

% Randomly choose trials for each bootstraps
for j = 1:options.j_for_decoder
    
    parfor_progress(options.bootstrapN)
    parfor nn = 1: options.bootstrapN
        % Random permute trials in every condition and every time windows
        % before training trials selection, so that each boot use different trials
%         pseudo_trial_pool_perm{j,nn} = cellfun(@(x) x(randperm(size(x,1)),1), pseudo_trial_decoder{j}, 'UniformOutput',false);
        pseudo_trial_pool_perm{j,nn} = cellfun(@(x) randperm(size(x,1)), pseudo_trial_decoder{j}(:,:,1), 'UniformOutput',false);   % Using same trials for the different temporal decoders in the same bootNum 
%         training_set{j,nn} = cellfun(@(x) x(1:options.min_reps4training,:),pseudo_trial_pool_perm{j,nn},'UniformOutput',false);
        training_set{j,nn} = cellfun(@(x) x(1:options.min_reps4training),pseudo_trial_pool_perm{j,nn},'UniformOutput',false);
        
        test_set_temp = cell(6,1);
        for c = 1: size(pseudo_trial_decoder{j},2)
%             test_set_temp{c} = cellfun(@(x) x(options.min_reps4training+1:min_reps_each_condition(c),:),pseudo_trial_pool_perm{j,nn}(:,c,:),'UniformOutput',false);
            test_set_temp{c} = cellfun(@(x) x(options.min_reps4training+1:min_reps_each_condition(c)),pseudo_trial_pool_perm{j,nn}(:,c,:),'UniformOutput',false);

        end
        
        testing_set{j,nn} = [test_set_temp{:}];
        
%         testing_set{j,nn} = [testing_set{j,nn}{:}];
        parfor_progress;
    end
end
parfor_progress(0);

%% Return
% decoder_prepare.FR = pseudo_trial_decoder;
% decoder_prepare.permutation_FR = pseudo_trial_pool_perm;
decoder_prepare.t_centers = decoder_tcenters;
decoder_prepare.cell_selected_ind = enough_reps_cell;
decoder_prepare.min_reps_each_condition =min_reps_each_condition;
decoder_prepare.min_reps4training = options.min_reps4training;
decoder_prepare.min_reps4each_condition = options.min_reps4each_condition;
decoder_prepare.decoder_window = options.decoder_window;
decoder_prepare.decoder_step_size = options.decoder_step_size;
decoder_prepare.bootstrapN = options.bootstrapN;
decoder_prepare.teacher_signal = answers_choice;
decoder_prepare.training_set = training_set;
decoder_prepare.testing_set = testing_set;
decoder_prepare.FR = pseudo_trial_decoder;  % added by ZZ @20231121

end
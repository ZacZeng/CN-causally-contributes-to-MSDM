%% Group microstimulation analysis from excel
% Begin from HH
% ZZ 20211015

% Read pre-analysed data from  Datahub and do some correlation analysis
% Positive here indicate contralateral preference (i.e. Rightward preference)

function function_handles = Group_MicroStim_xls(XlsData,if_tolerance)

% Read xls
num = XlsData.num;
txt = XlsData.txt;
raw  = XlsData.raw;
header = XlsData.header;

% Stim types
stim_type_list = {'Vestibular', 'Visual', 'Combined'};
stim_type_color = {'b', 'r', 'g'};
colorcode  = [41 89 204; 248 28 83; 14 153 46]/255;

% Comparison options
BIAS = 1; THRESHOLD = 2;
ttest_thresh = {0 1};  % t-test threshold

% == Figure default
set(0, 'defaultFigureRenderer', 'painters')
set(0, 'defaultFigureRendererMode', 'Manual')
set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46]/255);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======= Change the file directory===========

% Batch_address
mat_address = {
    % Major protocol goes here (Address, Suffix)
    'D:\Paper_rawdata\Raw_data\CN\Microstimulation\self-motion period\CN_m15&13_uStim', 'ustim'; 
    
    % Associative protocols
    'D:\Paper_rawdata\Raw_data\CN\Microstimulation\self-motion period\CN_m15&13_uStim_MemSac', 'MemSac';
    
    'D:\Paper_rawdata\Raw_data\CN\Microstimulation\self-motion period\CN_m15&13_uStim_Heading', 'PSTH';
    };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Global Mask for each protocol
mask_all = {
    strcmp(txt(:,header.Protocol), 'u-stim') & strcmp(txt(:,header.Area),'CD') & num(:,header.rep_ustim)>=8 & ...
    num(:,header.uA)<=80 & strcmp(txt(:,header.Note),'0:1500') & (num(:,header.Monkey) == 15 | num(:,header.Monkey) == 13 );
    
    
    strcmp(txt(:,header.Protocol),'MemSac') & strcmp(txt(:,header.Area),'CD')...
    & strcmp(txt(:,header.Unit),'MU21') & (num(:,header.Monkey) == 15 | num(:,header.Monkey) == 13 );
    
    
    strcmp(txt(:,header.Protocol),'HD') & strcmp(txt(:,header.Area),'CD')...
    & strcmp(txt(:,header.Unit),'MU21') & (num(:,header.Monkey) == 15 | num(:,header.Monkey) == 13 );
    
    };

% MonkeyID
monkey_name = {'Monkey D', 'Monkey F'};  % Duncan, Monkey
monkey_id = {15, 13};
monkey_shape = {'o', '^'};

monkey_included_for_loading = [15, 13];

monkey_mask_for_loading = false(size(num,1),1);
for mm = 1:length(monkey_included_for_loading)
    monkey_mask_for_loading = monkey_mask_for_loading | (num(:,header.Monkey) == monkey_included_for_loading(mm));
end

% Now apply monkey mask
for mm = 1:length(mask_all)
    mask_all{mm} = mask_all{mm} & monkey_mask_for_loading;
end

for pp = 1:size(mat_address, 1)
    xls_num{pp} = num(mask_all{pp},:);
    xls_txt{pp} = txt(mask_all{pp},:);
    xls_raw{pp} = raw(mask_all{pp}, :);
    
    % Basic information : E.g. '20140721m05s034h2x04y12d08478u6m5c67r2'
    cell_info_tmp = xls_num{pp}(:, [header.Date:header.Yloc header.Depth header.Chan1]);
    
    cell_info{pp} = strsplit(sprintf('%8dm%02ds%03dh%01dx%02dy%02dd%05du%02d\n', cell_info_tmp'),'\n')';
    cell_info{pp} = strcat(cell_info{pp}(1:end-1),xls_txt{pp}(:, header.FileNo));
    
end

%% Get Data
global group_Microstim_result;
% This lets me reload the whole m.file without group_result, which speeds up mu debugging
% Finnally added by ZZ @20230417

if isempty(group_Microstim_result)
    group_Microstim_result(size(xls_txt{1},1)).cellID = [];
    load_group_result();
end

    function load_group_result()
        
        tic;
        progressbar('Load .mat files');
        not_match = 0;
        
        %         for major_i = 1:length(cell_info{1})
        for major_i = 1:length(group_Microstim_result)    % Major protocol loop
            for pp = 1:size(mat_address,1)
                if pp ==1
                    % Get basic information for cellID
                    file_i = major_i;
                else
                    match_i = find(strncmp(cell_info{1}{major_i}, cell_info{pp}, 29));
                    
                    if length(match_i) == 1
                        file_i = match_i;
                    elseif length(match_i) > 1
                        fprintf('More than one match have been found for %s, \n Major ID = %s\n', xls_txt{1}{major_i, header.FileNo}, cell_info{1}{major_i});
                        file_i = match_i(1);  % The first appearance by default
                    else
                        if pp == 2
                            fprintf('No exact matched MemSac file for %s, ID = %s\n', xls_txt{1}{major_i,header.FileNo}, cell_info{1}{major_i});
                        elseif pp == 3
                            fprintf('No exact matched HeadingDiscrim file for %s, ID = %s\n', xls_txt{1}{major_i,header.FileNo}, cell_info{1}{major_i});
                        end
                        not_match = not_match +1;
                        file_i= NaN;
                    end
                end
                %             end
                
                % Load .mat for major and associative protocols
                if ~isnan(file_i)
                    try
                        mat_file_name = sprintf('%s_%g', xls_txt{pp}{file_i,header.FileNo},xls_num{pp}(file_i,header.Chan1));
                        mat_file_fullname = [mat_address{pp,1} '\' mat_file_name '_' mat_address{pp,2}];
                        
                        raw = load(mat_file_fullname);
                        
                        group_Microstim_result(major_i).cellID{pp} = cell_info{pp}{file_i};
                        group_Microstim_result(major_i).(['mat_raw_' mat_address{pp,2}]) = raw.result;
                        
                    catch
                        fprintf('Error Loading %s\n', [mat_file_name '_' mat_address{pp,2}]);
                        if pp == 1
                            keyboard
                        end
                    end
                    
                else
                    group_Microstim_result(major_i).(['mat_raw_' mat_address{pp,2}]) = [];
                end
                
            end
            
            progressbar(major_i / length(group_Microstim_result));
            
        end
        
        fprintf('Loaded %g files (%g of them are incomplete).\n', length(group_Microstim_result), not_match);
        toc;
    end


%% Construct more complete group_result for easier access
counts = 0;
for i = 1:length(group_Microstim_result)
    
    group_Microstim_result(i).hemisphere = xls_num{1}(i,header.Hemisphere);
    
    if isempty(group_Microstim_result(i).mat_raw_ustim)
        continue;
    end
    counts = counts+1;
    
    % 1) Microstim effect
    group_Microstim_result(i).ustim_repN = group_Microstim_result(i).mat_raw_ustim.repetitionN;
    group_Microstim_result(i).ustim_unique_heading = group_Microstim_result(i).mat_raw_ustim.unique_heading;
    group_Microstim_result(i).ustim_unique_condi = group_Microstim_result(i).mat_raw_ustim.unique_stim_type;
    group_Microstim_result(i).ustim_bias_psy = group_Microstim_result(i).mat_raw_ustim.Bias_psy;
    group_Microstim_result(i).ustim_p_bias = group_Microstim_result(i).mat_raw_ustim.P_bias;
    group_Microstim_result(i).ustim_bias_shift = group_Microstim_result(i).mat_raw_ustim.psy_bias_shift;
    group_Microstim_result(i).ustim_thresh_psy = group_Microstim_result(i).mat_raw_ustim.Thresh_psy;
    group_Microstim_result(i).ustim_p_thresh = group_Microstim_result(i).mat_raw_ustim.P_slope;
    group_Microstim_result(i).ustim_thresh_shift = group_Microstim_result(i).mat_raw_ustim.psy_thresh_shift;
    group_Microstim_result(i).ustim_correct_rate = group_Microstim_result(i).mat_raw_ustim.correct_rate;
    group_Microstim_result(i).ustim_correct_rate_shift = group_Microstim_result(i).mat_raw_ustim.correct_rate_shift;
    group_Microstim_result(i).bootstrap_bias_if_sig = group_Microstim_result(i).mat_raw_ustim.bootstrap_bias_if_sig;
    group_Microstim_result(i).bootstrap_p_bias = group_Microstim_result(i).mat_raw_ustim.bootstrap_p_bias;
    group_Microstim_result(i).bootstrap_thres_if_sig = group_Microstim_result(i).mat_raw_ustim.bootstrap_thres_if_sig;
    group_Microstim_result(i).bootstrap_p_thres = group_Microstim_result(i).mat_raw_ustim.bootstrap_p_thres;
    
    if ~isempty(group_Microstim_result(i).mat_raw_PSTH)
        % 2) PSTH from Heading task
        group_Microstim_result(i).psth_repN = group_Microstim_result(i).mat_raw_PSTH.repetitionN;
        group_Microstim_result(i).psth_unique_heading = group_Microstim_result(i).mat_raw_PSTH.CP{1,1}.Psy_func(:,1);
        group_Microstim_result(i).PREF_psth = group_Microstim_result(i).mat_raw_PSTH.PREF;
        group_Microstim_result(i).ChoicePreference = group_Microstim_result(i).mat_raw_PSTH.ChoicePreference;
        group_Microstim_result(i).ChoicePreference_pvalue = group_Microstim_result(i).mat_raw_PSTH.ChoicePreference_pvalue;
        
        % In TEMPO_GUI, choice preference uses the cell's PREF as its preferred
        % direction, because the hemisphere is unknow unless accessible to
        % Result.xls
        group_Microstim_result(i).if_pref_contralateral = group_Microstim_result(i).hemisphere ~= group_Microstim_result(i).PREF_psth;    % Defined by all correct trials
        group_Microstim_result(i).ChoicePreference = group_Microstim_result(i).ChoicePreference * sign(group_Microstim_result(i).if_pref_contralateral - 0.5);   % positive ChoicePref indicate contralateral preference
        
        % 3) Here comes the 'Modality preference'  (1-2, 1-3, 2-3)
        % May be i can compare whether microstimulation induces sites
        % preferring one modality shift more in this modal condition
        group_Microstim_result(i).ModalityPreference = group_Microstim_result(i).mat_raw_PSTH.ModalityPreference;
        group_Microstim_result(i).ModalityPreference_pvalue = group_Microstim_result(i).mat_raw_PSTH.ModalityPreference_pvalue;
    end
    
    % 4) MemSac dynamics
    if ~isempty(group_Microstim_result(i).mat_raw_MemSac)
        group_Microstim_result(i).MemSac_p = group_Microstim_result(i).mat_raw_MemSac.p;
        group_Microstim_result(i).MemSac_DDI = group_Microstim_result(i).mat_raw_MemSac.DDI;
        group_Microstim_result(i).MemSac_vectSum = aziToHeading(group_Microstim_result(i).mat_raw_MemSac.vectSum);  % negtive is left preferring
        if_contral = sign(group_Microstim_result(i).MemSac_vectSum) - group_Microstim_result(i).hemisphere;
        group_Microstim_result(i).MemSac_if_pref_contra = (if_contral==0 | if_contral==-3);
        group_Microstim_result(i).MemSac_AI = group_Microstim_result(i).mat_raw_MemSac.activityIndex;
    end
    
end

N = length(group_Microstim_result);
representative_cell = 40;  % Maybe I can let a site with significant shift in all conditions be here
unique_heading = group_Microstim_result(representative_cell).ustim_unique_heading;
%%  Initiation value
% Microstimulation
PSE_ctrl = NaN(N, 3);    Thresh_ctrl = NaN(N,3);
PSE_ustim = NaN(N, 3); Thresh_ustim = NaN(N,3);
dPSE = NaN(N, 3);          rThresh = NaN(N,3);      
pPSE = NaN(N,3);        pThresh = NaN(N,3);        
bootstrap_bias_if_sig = NaN(N,3);   bootstrap_thres_if_sig = NaN(N,3);  
bootstrap_p_bias = bootstrap_bias_if_sig;  bootstrap_p_thres = bootstrap_thres_if_sig; 
unique_heading_num = length(group_Microstim_result(representative_cell).ustim_unique_heading);
correct_rate_ctrl = NaN(N,unique_heading_num,3);  correct_rate_ustim = NaN(N,unique_heading_num,3);
% Shift, time window = 5 reps
shift_n = length(group_Microstim_result(representative_cell).ustim_bias_shift{1,1});
PSE_shift_ctrl = NaN(N,shift_n,3);  PSE_shift_ustim= NaN(N,shift_n,3);  % N sessions * 6 shift steps * 3 stim_type
thresh_shift_ctrl = NaN(N,shift_n,3); thresh_shift_ustim = NaN(N,shift_n,3);
correct_rate_shift_ctrl = NaN(N, unique_heading_num, shift_n, 3);    % N sessions * 9 headings * 6 shift steps * 3 stim_type
correct_rate_shift_ustim = NaN(N, unique_heading_num, shift_n, 3);
delta_psy_title ={'PSE Shift (o)'; 'Threshold Change (o)'};

% MemSac
memsac_DDI = NaN(N,6);  memsac_pref = NaN(N,6);
memsac_pvalue = NaN(N,6);
memsac_DDI_if_contra = NaN(N,6);

% Heading task
HD_ChDiv = NaN(N,3);
HD_CD_pvalue = NaN(N,3);
HD_ChDiv_if_contra = NaN(N,3);


%%  Microstimulation results from excel or Group_result
for i = 1:N
    
    % Microstimulation
    PSE_ctrl(i,:) = [group_Microstim_result(i).ustim_bias_psy{1,:}];
    PSE_ustim(i,:) = [group_Microstim_result(i).ustim_bias_psy{2,:}];
    pPSE(i,:) = group_Microstim_result(i).ustim_p_bias;
    bootstrap_bias_if_sig(i,:) = group_Microstim_result(i).bootstrap_bias_if_sig;
    bootstrap_p_bias(i,:) = group_Microstim_result(i).bootstrap_p_bias;
    
    Thresh_ctrl(i,:) = [group_Microstim_result(i).ustim_thresh_psy{1,:}];
    Thresh_ustim(i,:) = [group_Microstim_result(i).ustim_thresh_psy{2,:}];
    pThresh(i,:) = group_Microstim_result(i).ustim_p_thresh;
    bootstrap_thres_if_sig(i,:) = group_Microstim_result(i).bootstrap_thres_if_sig;
    bootstrap_p_thres(i,:) = group_Microstim_result(i).bootstrap_p_thres;
    
    correct_rate_ctrl(i,:,:) = (group_Microstim_result(i).ustim_correct_rate{1})';       % N sessions * 9 headings * 3 stim_type
    correct_rate_ustim(i,:,:) = (group_Microstim_result(i).ustim_correct_rate{2})';
    
    % Temporal Dynamic changes with trials going
    if group_Microstim_result(i).ustim_repN == 10      % For simplicity, only choose 10 reps blocks
        PSE_shift_ctrl(i,:,:) = cell2mat([group_Microstim_result(i).ustim_bias_shift(:,1)])';  % N sessions * 6 shift steps * 3 stim_type
        PSE_shift_ustim(i,:,:) = cell2mat([group_Microstim_result(i).ustim_bias_shift(:,2)])';
        thresh_shift_ctrl(i,:,:) = cell2mat([group_Microstim_result(i).ustim_thresh_shift(:,1)])';
        thresh_shift_ustim(i,:,:) = cell2mat([group_Microstim_result(i).ustim_thresh_shift(:,2)])';
        correct_rate_shift_ctrl(i,:,:,:) = reshape([group_Microstim_result(i).ustim_correct_rate_shift{:,1}], unique_heading_num,[], 3);   % N sessions * 9 headings * 6 shift steps * 3 stim_type
        correct_rate_shift_ustim(i,:,:,:) = reshape([group_Microstim_result(i).ustim_correct_rate_shift{:,2}], unique_heading_num,[], 3);   % N sessions * 9 headings * 6 shift steps * 3 stim_type
    end
    
    % Switch microstimulation result according to hemisphere
    hemisphere(i) = group_Microstim_result(i).hemisphere;
    if hemisphere(i) == 2
        % Unify the Y-axis of psychometric curves as Contralateral choice and X-axis as from ipsi to conta headings
        PSE_ctrl(i,:) = -PSE_ctrl(i,:);
        PSE_ustim(i,:) = -PSE_ustim(i,:);
        PSE_shift_ctrl(i,:,:) = -PSE_shift_ctrl(i,:,:);
        PSE_shift_ustim(i,:,:) = -PSE_shift_ustim(i,:,:);
        
        % From the left to right is correct rate of ipsi to contra headings
        correct_rate_ctrl(i,:,:) = fliplr(correct_rate_ctrl(i,:,:));
        correct_rate_ustim(i,:,:) = fliplr(correct_rate_ustim(i,:,:));
    end
    
    % MemSac
    if ~isempty(group_Microstim_result(i).MemSac_DDI)
        memsac_DDI(i,:) = group_Microstim_result(i).MemSac_DDI;
        memsac_pvalue(i,:) = group_Microstim_result(i).MemSac_p;
        memsac_pref(i,:) = group_Microstim_result(i).MemSac_vectSum;    % negtive is left preferring
        memsac_DDI_if_contra(i,:) = group_Microstim_result(i).MemSac_if_pref_contra;
    end
    
    % Heading task
    j = 1;  % only want stimulus period
    if ~isempty(group_Microstim_result(i).ChoicePreference)
        HD_ChDiv(i,:) = group_Microstim_result(i).ChoicePreference(j,:);
        HD_CD_pvalue(i,:) = group_Microstim_result(i).ChoicePreference_pvalue(j,:);
        HD_ChDiv_if_contra(i,:) = group_Microstim_result(i).if_pref_contralateral;
    end
    
end

dPSE = PSE_ctrl - PSE_ustim;   % negtive, ipsilateral shift; positive, contralater shift
% % Switch sign of dPSE when right hemisphere to let
% % negative value indicates ipsilateral shift and positive indicate
% % contralater shift
% dPSE = -(dPSE .* sign(hemisphere - 1.5)');    % negtive, ipsilateral shift; positive, contralater shift
rThresh = -(Thresh_ustim - Thresh_ctrl);            % larger than 0 indicate worse sensitivity
dPSE_dynamic = PSE_shift_ctrl - PSE_shift_ustim;

psycho_ctrl = [PSE_ctrl Thresh_ctrl];
psycho_ustim = [PSE_ustim Thresh_ustim];
delta_psycho = [dPSE  rThresh];
p_value_psycho = [pPSE pThresh];
bootstrap_if_sig = [bootstrap_bias_if_sig bootstrap_thres_if_sig]; 
bootstrap_p_value = [bootstrap_p_bias bootstrap_p_thres]; 
%% Cell Selection and Cell Counter
% Now we can get data from xls. as well as group_result.mat
cell_selection();
    function cell_selection(t_cell_selection_num)
        
        if nargin < 1
            t_cell_selection_num = 2;
        end
        
        % Limitations on repetition number
        select_all_all_monkey = ([group_Microstim_result.ustim_repN]' >= 8);
        
        select_bottom_line_all_monkey = select_all_all_monkey;
        
        % + T(ypical) Sites
        
        t_cell_selection_criteria = ...
            {  %  Logic                                                Notes
            'Bottom-line (all )', ones(size(select_all_all_monkey));
            % u_stim based
            'ustim p value (any)', any(bootstrap_p_bias<0.01 , 2);
             'ustim threshold p value (any)', any(bootstrap_p_thres<0.01 , 2);
            % Mem-Sac based
            'Memory p value', memsac_pvalue(:,3) < 0.05;
            % HD based
            'ChoicePreference p value (any)', any(HD_CD_pvalue<0.01, 2);
            % ustim and HD 
            % Added @ 20230726 refer to Ding 2012 neuron Fig 4
            
            % Hemisphere based
            'Left Hemisphere', hemisphere' == 1;   % almost all sites are from left hemisphere
            };
        
        select_tsites_all_monkey = select_bottom_line_all_monkey & t_cell_selection_criteria{t_cell_selection_num, 2};
        select_no_tsites_all_monkey = select_bottom_line_all_monkey & ~t_cell_selection_criteria{t_cell_selection_num, 2};
        t_criterion_txt = t_cell_selection_criteria{t_cell_selection_num, 1};
        
        % ------- Count cell numbers for each monkey -----------------
        n_monkey = length(monkey_included_for_loading);
        cell_nums = zeros(n_monkey+1, 3);   % All sites / significant sites / non-sig sites
        for mm = 1:n_monkey
            select_monkey{mm} = xls_num{1}(:, header.Monkey) == monkey_included_for_loading(mm);
            cell_nums(mm,1) = sum(select_all_all_monkey & select_monkey{mm});
            cell_nums(mm,2) = sum(select_tsites_all_monkey & select_monkey{mm});
            cell_nums(mm,3) = sum(select_no_tsites_all_monkey & select_monkey{mm});
        end
        
        % ------------- Update actual dataset for analysis ------------
        monkey_included_for_analysis = monkey_included_for_loading(logical([get(findall(gcbf,'tag','Duncan_data'),'value') get(findall(gcbf,'tag','Fara_data'),'value')]));
        monkey_mask_for_analysis = false(length(group_Microstim_result),1);
        for mm = 1:length(monkey_included_for_analysis)
            monkey_mask_for_analysis = monkey_mask_for_analysis | (xls_num{1}(:, header.Monkey) == monkey_included_for_analysis(mm));
        end
        
        % -----------------  Affect all analysis below --------------------
        select_all = select_all_all_monkey & monkey_mask_for_analysis;
        select_tsites = select_tsites_all_monkey & monkey_mask_for_analysis;
        select_no_tsites = select_no_tsites_all_monkey & monkey_mask_for_analysis;
        
        cell_nums(end,:) = [sum(select_all) sum(select_tsites) sum(select_no_tsites)];
        
        % ------------------------Update Cell Counter ------------------
        h_all = findall(gcbf,'tag','num_all_units');
        set(h_all, 'string', sprintf('%7d%7d%7d\n',cell_nums'), 'fontsize',13);
        h_t_criterion = findall(gcbf,'tag','t_criterion');
        set(h_t_criterion, 'string',{t_cell_selection_criteria{:,1}});
        set(h_t_criterion, 'value',t_cell_selection_num);
    end

%%
    function cell_position()
        
        GM = -1;
        CD = -3;
        CC = -4;
        NAc = -5;
        GP = -6;
        Pu = -7;
        FC = -8;
        LV = -9;
        
        drawmapping_data = { % [Monkey,hemi]    % Grid No. of (AP0,Middle)   % Data
            
        [15,1],[31,0]     % Duncan left
        [15,2],[31,0]
        [13,1],[19,0]     % Fara left
        [13,2],[20,0]
        };
    
    drawmapping_data{1,3} = {
        %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] , electrode retrieval}
        
        {127, [28, 10],[2.3,0.0], [CD 62 95]}       %%%%%%  Microstimulation From here
        {129, [32, 9],[2.0 0.0], [GM 0 16; CD 93 115]}
        {130, [30,8], [1.8 0.0], [GM 11 20; CD 119 132]}
        {131, [17,9], [1.8,0.0], [GM 2 23; CD 129 150]}
        {132, [20,10], [2.0,0.0],[CD 101 125]} % only to 12243
        {133, [18,10],[2.0,0.0],[GM 0 22; CD 105 118]} %only to11719
        {134, [19,10],[2.0,0.0], [CD 100 120]} % only to 11548
        {135, [34,10],[2.0,0.0], [GM 2 35; CD 94 110]} % only to 10948
        {136, [34,10],[2.0,0.0], [CD 94 110]}  % only to 10780
        {137, [35,10],[1.8,0.0],[CD 128 142]} % only to 14187
        {138, [30,10],[2.0,0.0],[GM 14 21; CD 106 124]} % only to 12366
        {139, [25,9],[2.0,0.0],[CD 102 130]}
        {140, [24,10],[2.0,0.0],[GM 0 24; CD 99 120]} % only to 11942
        {141, [26,10],[2.0,0.0],[GM 6 15; CD 104 121]} % only to 12067
        {142, [29,10],[2.0,0.0],[CD 97 121]}
        {143, [29, 8],[2.0,0.0],[CC 25 50; CD 86 100]} % only to 9723
        {144, [27, 8],[2.0,0.0],[CC 33 56; CD 97 110]} % only to 10882
        {145, [28,7], [1.8, 0.0], [CC 58 81;CD 122 135]} % only to 13346
        {146, [25, 8],[2.0,0.0], [CD 104 132]}
        {147, [28,7],[2.0,0.0],[CC 27 54; CD 94 105]} % only to 10352
        {148, [27,9],[2.0,0.0],[CD 91 110]} % only to 10763
        {149, [26,10],[2.0,0.0],[CD 110 138]}
        {150, [24,8],[1.8,0.0],[GM 0 23; CD 110 135]} % only to 13420
        {151, [25,7],[1.8,0.0],[GM 2 26; CD 117 135]} % only to 13384
        {152, [30,9],[1.8,0.0],[GM 0 23; CD 111 125]} % only to 12140
        {153, [22,10],[2.0,0.0],[GM 0 18; GM 42 60; CD 110 120]} % only to11488
        {154, [20,9],[2.0,0.0],[CD 98 120]} % only to 11760
        {155, [18,8],[1.9,0.0],[GM 1 19; CD 119 146 ]}
        {156, [17,10],[2.0,0.0],[CD 122 143]}
        {157, [30,7],[1.8,0.0],[CC 55 70; CD 124 150]}
        {158, [32,10],[1.9,0.0],[GM 0 19; CD 109 138]}
        {159, [31,10],[1.8,0.0],[CD 122 145]}
        {160, [26,8],[1.8,0.0],[CD 132 150]}
        {161, [22,9],[2.0,0.0],[CD 106 120]} %only to 11878
        {162, [28,9],[1.9,0.0],[CD 102 125]} % only to 12049
        {164, [28,9],[2.0,0.0],[CD 93 120]} % only to 11171
        {165, [31,9],[1.9,0.0],[CD 111 125]} % only to 125
        {166, [31,9],[1.9,0.0],[CD 111 135]}
        {167, [27,11],[1.9,0,0],[CD 110 130]}
        {168, [24,7],[1.9,0.0],[GM 6 25; CD 114 136]}
        {169, [31,9],[1.9,0.0],[GM 3 25; CD 111 140]}
        {170, [32,8],[1.8,0.0],[GM 0 16; CC 57 73; CD 115 145]}
        {171, [29,9],[1.9,0.0],[GM 5 23; CD 108 125]} % only to 12329
        {172, [26,9],[1.9,0.0],[GM 0 19; CD 105 120]} % only to 11510
        {173, [23,10],[2.0,0.0],[GM 0 18; GM 48 62; CD 100 135]}
        {174, [20,8],[2.0,0.0],[CD 106 130]}
        {175, [24,9],[2.0,0.0],[CD 104 125]}
        {176, [25,10],[2.0,0.0],[CD 102 110]} % only to 10550
        {177, [25,10],[2.0,0.0],[CD 102 120]} % only to 11553
        {178, [27,10],[2.0,0.0],[CD 101 110]}
        {179, [27,10],[1.9,0.0],[GM 7 28; CD 113 125 ]} % only to 12249
        {180, [27,10],[1.9,0.1],[CD 101 120]}
        {181, [31,8],[1.9,0.0],[CD 120 125]} % only to 12473
        {182, [31,8],[1.9,0.0],[CC 44 71; CD 120 150]}
        {183, [33,10],[1.9,0.0],[GM 2 25; CD 110 125]} % only to 12085
        {184, [33,9],[1.9,0.0],[CD 113 120]} % only to 11803
        {185, [33,9],[1.9,0.0],[GM 7 17; CD 113 135]}
        {186, [34,9],[1.9,0.0],[GM 0 22; CD 114 125]} % only to 12397
        {187, [31,7],[1.8,0.0],[CC 58 87; CD 128 150]}
        {188, [33,7],[1.9,0.0],[GM 6 30; CC 43 61; CC 65 82; CD 121 150]}
        {189, [33,8],[1.9,0.0],[GM 3 15; CC 55 73; CD 115 140]}
        {190, [28,8],[1.9,0.0],[GM 7 21; CD 116 140]}
        {191, [21,9],[1.9,0.0],[GM 7 26; CD 115 130]}; %only to 12645
        {192, [21,10],[1.9,0.0],[GM 8 35; CD 122 140]}
        
        }';
    
    
    drawmapping_data{3,3} = {
        {272, [15, 8], [2.3, 0.0], [GM 8 35; CD 141  165]}  % only to 165
        {273, [15, 6], [2.3, 0.0], [GM 0 26; GM 68 97; CD 142 170]} % only to 165
        {274, [15, 10], [2.4, 0.0], [GM 13 35; CD 132 145]}   % only to 143
        {275, [15, 9], [2.4, 0.0], [GM 8 35; CD 134 150]}    % only to 150
        {276, [15, 7], [2.3, 0.0], [GM 15 35; CD 137 160]} % only to 160
        {278, [17, 8], [2.3, 0.0], [GM 10 30; CD 145 160]} % only to 160
        {279, [19, 8], [2.3, 0.0], [GM 8 30; CD 142 180]} % only to 175
        {280, [16, 8], [2.3, 0.0], [GM 11 39; CD 147 180]} % only to 176
        {281, [16, 7], [2.4, 0.0], [GM 14 32; CD 138 170]} % only to 163
        {282, [16, 9], [2.4, 0.0], [CD 145 170]} % only to 165
        {283, [16, 6], [2.3, 0.0], [GM 12 36; GM 78 97; CD 154 185]} % only to 181
        {284, [17, 9], [2.3, 0.0], [GM 22 43; CD 159 195]} % only to 190
        {285, [18, 7], [2.3, 0.0], [GM 20 39; CD 158 175]} % only to 173
        {286, [18, 8], [2.3, 0.0], [GM 80 92; CD 150 165]} % only to 160
        {287, [18, 6], [2.3, 0.0], [GM 67 95; CD 146 170]} % only to 169
        {288, [18, 9], [2.3, 0.0], [GM 12 38; CD 145 165]} % only to 162
        {289, [19, 9], [2.3, 0.0], [GM 13 31; GM 67 100; CD 148 160]} % only to 158
        {290, [19, 7], [2.4, 0.0], [GM 0 20; GM  57 83; CD 138 165]} % only to 160
        {291, [20, 7], [2.3, 0.0], [GM 15 31; CD 148 170]} % only to 168
        {292, [20, 8], [2.3, 0.0], [GM 13 42; CD 145 165]} % only to 160
        {293, [20, 6], [2.3, 0.0], [GM 7 20; GM 64 96; CD 150 182; GM 211 217]}
        {294, [21, 8], [2.3, 0.0], [GM 13 35; GM 78 93; CD 156 180]} % only to 175
        {295, [21, 7], [2.3, 0.0], [GM 11 27; GM 66 102; CD 159 186]}
        {296, [21, 9], [2.3, 0.0], [GM 16 39; CD 146 165]}   % only to 160
        {297, [21, 10], [2.3, 0.0], [GM 19 47; CD 151 170]}
        {298, [19, 10], [2.3, 0.0], [GM 16 39; CD 146 175]}
        {299, [22, 7], [2.3, 0.0], [GM 23 31; GM 73 107; CD 157 190]}
        {300, [22, 8], [2.3, 0.0], [GM 15 33; CD 151 190]}
        {301, [17, 7], [2.3, 0.0], [GM 81 97; CD 156 170]}  % only to 168
        {302, [17, 10], [2.4, 0.0], [GM 81 97; CD 155 171]}
        {303, [17, 6], [2.4, 0.0], [GM 81 96; CD 148 180]}  % only to 176
        {304, [14, 8], [2.4, 0.0], [GM 33 52; CD 156 175]}  % only to 172
        {305, [14, 7], [2.4, 0.0], [GM 17 45; CD 154 175]}  % only to 173
        {306, [14, 6], [2.4, 0.0], [GM 27 49; CD 152 160]}  % only to 156
        {307, [14, 9], [2.4, 0.0], [GM 16 32; CD 146 155]}  % only to 152
        {308, [13, 8], [2.4, 0.0], [GM 24 40; CD 150 165]}  % only to 163
        {309, [13, 9], [2.4, 0.0], [GM 22 41; CD 152 165]}  % only to 161
        {310, [13, 7], [2.4, 0.0], [GM 18 38; CD 147 160]}  % only to 159
        {311, [13, 6], [2.4, 0.0], [GM 11 36; CD 147 165]}  % only to 161
        {312, [12, 8], [2.4, 0.0], [GM 12 36; CD 149 165]}  % only to 160
        {313, [12, 6], [2.4, 0.0], [GM 17 42; CD 146 150]}  % only to 148
        {314, [10, 7], [2.4, 0.0], [GM 12 37; CD 149 170]}  % only to 168
        {315, [11, 7], [2.4, 0.0], [GM 13 41; CD 148 160]}  % only to 154
        {316, [12, 7], [2.4, 0.0], [GM 15 40; CD 153 170]}  % only to 169
        {317, [12, 5], [2.4, 0.0], [GM 12 47; CD 156 175]}  % only to 172
        {318, [11, 8], [2.4, 0.0], [GM 25 48; CD 157 180]}  % only to 177
        {319, [11, 6], [2.4, 0.0], [GM 14 42; CD 146 155]}  % only to 153
        {320, [11, 9], [2.4, 0.0], [GM 17 41; GM 79 92; CD 151 160]}  % only to 156
        {321, [10, 6], [2.4, 0.0], [GM 17 37; CD 149 180]}  % only to 175
        {322, [10, 8], [2.4, 0.0], [GM 22 37; CD 154 175]}  % only to 172
        {323, [10, 5], [2.4, 0.0], [GM 20 41; GM 75 102; CD 151 180]} % only to 175
        {324, [9, 8], [2.4, 0.0], [GM 22 39; GM 89 107; CD 155 180]}  % only to 178
        }';
    
    
    drawmapping_data{4,3} = {
        {258, [10, 6], [1.9, 0.0], [GM 40 57; GM 123 143; GM 180 193]}  % only to 193
        {261, [15, 7], [2.3, 0.0], [GM 6 18; GM 71 80; CD 133 190]}
        {265, [15, 9], [2.3, 0.0], [GM 0 23; CD 130 166]}
        {267, [15, 10], [2.3, 0.0], [GM 8 26; GM 69 87; CD 133 167]}
        {268, [17, 7], [2.3, 0.0], [GM 8 27; CD 144 172]} % only to 172
        {269, [10, 8], [2.3, 0.0], [GM 6 24; CD 135 172]} % only to 171
        {270, [13, 8], [2.3, 0.0], [GM 5 12; CD 136 165]} % only to 164
        {271, [19, 7], [2.3, 0.0], [GM 57 87; CD 145 165]} % only to 163
        };
    
    for i = 1:length(cell_info{1})
        % Decode cell position
        thisID = sscanf(cell_info{1}{i}, '%gm%gs%gh%gx%gy%gd%g');
        this_monkey_hemi = thisID([2 4])';
        this_session = thisID(3);
        this_pos_raw = thisID(5:7)';
        
        % Find entry in drawmapping_data
        found = 0;
        for mmhh = 1:size(drawmapping_data,1)
            if all(drawmapping_data{mmhh,1} == this_monkey_hemi)
                xyoffsets = drawmapping_data{mmhh,2};
                this_mapping_data = drawmapping_data{mmhh,3};
                
                for ee = 1:length(this_mapping_data)
                    
                    if all(this_mapping_data{ee}{1}==this_session) && all(this_mapping_data{ee}{2} == this_pos_raw(1:2))
                        found = found +1;
                        
                        AP = (this_pos_raw(1) - xyoffsets(1)) *0.8;
                        
                        which_GM_is_CD = find(this_mapping_data{ee}{4}(:,1) == CD);
                        which_CD_this_cell_in = find((this_mapping_data{ee}{4}(which_GM_is_CD,2) <= ceil(this_pos_raw(3)/100)) & ...
                            (this_mapping_data{ee}{4}(which_GM_is_CD,3) >= fix(this_pos_raw(3)/100)));
                        
                        if which_CD_this_cell_in == 1   % should always equal 1, for there is no gap inside CD
                            this_surface = this_mapping_data{ee}{4}(which_GM_is_CD(which_CD_this_cell_in), 2);
                            depth = this_pos_raw(3) - 100 * this_surface;
                        elseif which_CD_this_cell_in == 2
                            depth = nan;
                        else
                            disp('Check CD Layers !!!'); beep; keyboard;
                        end
                        
                        site_position(i,:) = [this_monkey_hemi AP depth which_CD_this_cell_in] ;
                        
                    end
                end
            end
            
        end
        
    end
    
    end

%% ====================== Function Handles ====================================
function_handles = [];
function_handles = {
    'Microstimulation Effects', {
    'Bias and Threshold comparison across sessions', @f1p1;
    'Bias and Threshold change (Bar)', @f1p2;
    'Bias and Threshold change comparison across conditions', @f1p3;
    'Normalized Bias comparison across conditions', @f1p4;
    'Correct rate of each heading', @f1p5;
    'Temporal dynamics of dPSE', @f1p6;
    };
    
    'Correlations', {
    'HeatMap in Three task', @f2p1
    'MemSac with dPSE' , @f2p2
    'ChoiceDive with dPSE' , @f2p3
    'MemSac with ChoiceDive', @f2p4
    };
    
    'Site Distribution', {
    'Microstim-site Distribution' , @f3p1
    'Site-distribution across monkeys', @f3p2
    };
    
    'NoShow',{@cell_selection};
    
    };
% ==========================================================================
%%
% Example initiation for site selection
methods_of_select = {
    select_all , 'All sites'
    select_tsites , 'Typical sites'
    select_no_tsites, 'Non-typical sites'
    };
monkey_included_for_analysis;
select_monkey;

figN = 1021;

%%
    function f1p1(debug)       % Session-by-session bias and threshold comparison
        if debug; dbstack; keyboard; end
        %% Microstimulation effects in all conditions and comparison between them
        % Monkey mask has been included
        methods_of_select = {
            select_all , 'All sites'
            select_tsites , 'Typical sites'
            select_no_tsites, 'Non-typical sites'
            };
        
        %         %  Scripts written before, called here for simplicity
        %         Group_Microstim_ZZ();    % PSE and threshold comparison, no DDI
        
        for ms = 1: size(methods_of_select, 1)
            
            % Median value
            psycho_ctrl_median{ms} = nanmedian(psycho_ctrl(methods_of_select{ms,1},:));
            psycho_ustim_median{ms} = nanmedian(psycho_ustim(methods_of_select{ms,1},:));
            delta_psycho_median{ms} = nanmedian(delta_psycho(methods_of_select{ms,1},:));
            
            % Sorted by method of selection
            psycho_ctrl_select{ms} = psycho_ctrl(methods_of_select{ms,1},:);
            psycho_ustim_select{ms} = psycho_ustim(methods_of_select{ms,1},:);
            p_value_psycho_select{ms} = p_value_psycho(methods_of_select{ms,1},:);
            
            % == Plotting
            set(figure(figN), 'name', ['Session-by-session comparison' methods_of_select{ms, 2}], 'pos', [20 10 1000 1000]); clf; figN=figN+1;
            hb = tight_subplot(3,2,[.1 .05],[.1 .1], .1);
            
            for option = BIAS : THRESHOLD
                for condi = 1:3
                    p = signrank(psycho_ctrl_select{ms}(:, condi+(option-1)*3), psycho_ustim_select{ms}(:, condi+(option-1)*3));     % Wilcoxon signed rank test
                    
                    set(gcf, 'CurrentAxes', hb(condi+(option-1)*3)); hold on; 
                    
                    % For marking dots with significant effects
%                     sig_ind = p_value_psycho(:, condi+(option-1)*3) <= 0.05;
%                     sig_ind = bootstrap_if_sig(:, condi+(option-1)*3);
                    sig_ind = bootstrap_p_value(:, condi+(option-1)*3) <= 0.01; 
                    
                    for mon = monkey_included_for_analysis
                        m = find(monkey_included_for_loading==mon);
                        
                        % == Plotting
                        plot(psycho_ctrl(methods_of_select{ms,1}&select_monkey{m}, condi+(option-1)*3),...
                            psycho_ustim(methods_of_select{ms,1}&select_monkey{m}, condi+(option-1)*3),...
                            monkey_shape{m}, 'Color', colorcode(condi,:), 'MarkerSize',8,'LineWidth',0.75);
                        plot(psycho_ctrl(methods_of_select{ms,1}&sig_ind&select_monkey{m}, condi+(option-1)*3),...
                            psycho_ustim(methods_of_select{ms,1}&sig_ind&select_monkey{m}, condi+(option-1)*3),...
                            monkey_shape{m}, 'Color', colorcode(condi,:), 'MarkerSize',8, 'MarkerFaceColor', colorcode(condi,:),'LineWidth',0.75);
                    end
                    
                    axis square;

                    % Line of unity
                    min_scale = min([psycho_ctrl_select{ms}(:, condi+(option-1)*3);psycho_ustim_select{ms}(:,condi+(option-1)*3)]);
                    max_scale = max([psycho_ctrl_select{ms}(:, condi+(option-1)*3);psycho_ustim_select{ms}(:,condi+(option-1)*3)]);
                    axis([min_scale-0.1 max_scale+0.1 min_scale-0.1 max_scale+0.1]);
                    plot([min_scale-0.2 max_scale*1.05], [min_scale-0.2 max_scale*1.05], 'k--');
%                     plot([0 0], ylim, 'k--'); plot(xlim, [0 0], 'k--');
                    
                    
                    % Annotate median value
                    plot(psycho_ctrl_median{ms}(:, condi+(option-1)*3),min_scale+0.1,'v',...
                        'Color', colorcode(condi,:));
                    text(psycho_ctrl_median{ms}(:, condi+(option-1)*3),min_scale+0.1,...
                        num2str(psycho_ctrl_median{ms}(:, condi+(option-1)*3)),...
                        'Color', colorcode(condi,:), 'FontSize', 15, 'FontWeight', 'bold')     % Text location is always changed after axis square !!

                    plot(min_scale+0.1, psycho_ustim_median{ms}(:, condi+(option-1)*3),'<',...
                        'Color', colorcode(condi,:));
                    text(min_scale+0.1, psycho_ustim_median{ms}(:, condi+(option-1)*3),...
                        num2str(psycho_ustim_median{ms}(:, condi+(option-1)*3)),...
                        'Color', colorcode(condi,:), 'FontSize', 15, 'FontWeight', 'bold')
                    
                    % Annotate  N number and p-value
                    txt = {[ 'n = ', num2str(sum(~isnan(psycho_ctrl_select{ms}(:,condi+(option-1)*3)))) ], ['p: ', num2str(roundn(p,-4))], ['n_{sig} = ', num2str(sum(sig_ind&select_monkey{m})) ]};
                    str = text(min_scale-0.2+0.05,max_scale, txt, 'HorizontalAlignment','left');
                    set(str, 'Color', [p<0.05,(p<0.05)*0.7,0]);
                                        
                    % Title
                    if condi == 1
                        if option ==BIAS
                            title('Bias','Color', 'k', 'FontSize',20,'FontWeight','bold');
                        else
                            title('Threshold','Color','k', 'FontSize',20,'FontWeight','bold');
                        end
                    end
                    
                    if condi ==2
                        ylabel('\muStim (o)');
                    end
                    
                    if condi == 3
                        xlabel('No \muStim (o)');
                    end
                                        
                end
            end
            SetFigure();
        end
        
        
    end

%%
    function f1p2(debug)       % Bias and threshold change in each conditions (Bar plot)
        if debug; dbstack; keyboard; end
        
        % Monkey mask has been included
        methods_of_select = {
            select_all , 'All sites'
            select_tsites , 'Typical sites'
            select_no_tsites, 'Non-typical sites'
            };
        
        for ms = 1: size(methods_of_select, 1)
            % Median value
            delta_psycho_median{ms} = nanmedian(delta_psycho(methods_of_select{ms,1},:));
            
            % Sorted by method of selection
            delta_psycho_select{ms} = delta_psycho(methods_of_select{ms,1},:);
%             p_value_psycho_select{ms} = p_value_psycho(methods_of_select{ms,1},:);
            p_value_psycho_select{ms} = bootstrap_p_value(methods_of_select{ms,1},:);

        end
        
        % == Plotting
        nbin = 10;  nticks = 5;
        for ms = 1:size(methods_of_select, 1)
            set(figure(figN), 'name', ['CD_Microstim_Population' methods_of_select{ms,2}], 'Position', [40 10 1000 1000]); clf; figN=figN+1;
            hd = tight_subplot(3,2,[.1 .05],[.1 .1], .1);
            
            for option = BIAS : THRESHOLD
                for condi = 1: 3
                    % Whether significant at populational level
                    p = signrank(delta_psycho_select{ms}(:,condi+(option-1)*3));    % Wilcoxon signed rank test
                    
                    % Significant indication
                    Sig = p_value_psycho_select{ms}(:,condi+(option-1)*3) <= 0.01;
                    NSig = p_value_psycho_select{ms}(:,condi+(option-1)*3) > 0.01;
                    
                    axes(hd(condi+(option-1)*3));  hold on;
                    
                    [~, xcenters]  = hist(delta_psycho_select{ms}(:, condi+(option-1)*3), nbin);
                    % Centralizing
                    max_abs_xcenters = max(abs(xcenters));
                    xcenters = linspace(-max_abs_xcenters, max_abs_xcenters, nbin);
                    
                    histS= hist(delta_psycho_select{ms}(Sig, condi+(option-1)*3), xcenters);
                    histNS = hist(delta_psycho_select{ms}(NSig,condi+(option-1)*3), xcenters);
                    
                    hbars = bar(gca, xcenters, [histS' histNS'],1,'stacked','LineWidth',2);
                    set(hbars,'EdgeColor',colorcode(condi,:),'FaceColor',colorcode(condi,:));
                    set(hbars(2),'FaceColor','none');
                    
                    
                    % Annotation
                    if condi==3 && option==1
                        xlabel('PSE Shift (o)', 'FontWeight', 'bold');
                    end
                    
                    if condi==3 && option==2
                        xlabel('Threshold Change (o)','FontWeight', 'bold');
                    end
                    
                    if condi==2 && option==1
                        ylabel('Number of case','FontWeight', 'bold');
                    end
                    
                    xlim([min(xcenters), max(xcenters)]*1.2);
                    ylim([0, max(histS+histNS)+1]);
                    xticks(linspace(roundn(min(xcenters),-1),roundn(max(xcenters),-1),nticks));
                    
                    if option == 2
                        txt =['n = ', num2str(sum(~isnan(delta_psycho_select{ms}(:,condi))))];
                        text(max(xcenters), max(histS+histNS), txt, 'Color', colorcode(condi,:), 'HorizontalAlignment', 'right');
                    end
                    
                    plot(delta_psycho_median{ms}(condi+(option-1)*3)*[1 1], [0 max(histS+histNS)+1],...
                        'Color', [p<0.05,0.7*(p<0.05),0],'LineWidth', 2);
                    str = text(xcenters(3), max(histS+histNS)+1.5,...
                        sprintf('Median: %2.2f   p: %0.3f', delta_psycho_median{ms}(condi+(option-1)*3), p));
                    set(str, 'Color', [p<0.05,0.7*(p<0.05),0]);
                    
                end
            end
            
        end
    end

%%
    function f1p3(debug)       % Comparing Bias and threshold change between different conditions
        if debug; dbstack; keyboard; end
        
        % Monkey mask has been included
        methods_of_select = {
            select_all , 'All sites'
            select_tsites , 'Typical sites'
            select_no_tsites, 'Non-typical sites'
            };
        
        comp_pair = {[1 2]; [1 3]; [2 3]};   % comparing pairs
        
        % Significance indication
%         sig_ind = p_value_psycho <= 0.05;
        sig_ind = bootstrap_p_value <= 0.01;
        
        for ms = 1:size(methods_of_select, 1)
            % Median value
            delta_psycho_median{ms} = nanmedian(delta_psycho(methods_of_select{ms,1},:));
            % Selected value
            delta_psycho_select{ms} = delta_psycho(methods_of_select{ms,1},:);
            p_value_psycho_select{ms} = p_value_psycho(methods_of_select{ms,1},:);
            
            % -- Plotting
            set(figure(figN), 'name', ['Across-condition comparison ' methods_of_select{ms, 2}], 'Position',[80 10 1000 1000]); clf; figN=figN+1;
            ha = tight_subplot(3,2,[.1 .05],[.1 .1], .1);
            
            for option = BIAS:THRESHOLD
                for cp = 1:length(comp_pair)
                    axes(ha(cp+(option-1)*3)); hold on;
                    
                    Sig1 = sig_ind(:, comp_pair{cp}(1)+(option-1)*3) & ~sig_ind(:, comp_pair{cp}(2)+(option-1)*3);   % Only significant in the first condition
                    Sig2 = sig_ind(:, comp_pair{cp}(2)+(option-1)*3) & ~sig_ind(:, comp_pair{cp}(1)+(option-1)*3);   % Only significant in the second condition
                    Sig_both = sig_ind(:, comp_pair{cp}(2)+(option-1)*3) & sig_ind(:, comp_pair{cp}(1)+(option-1)*3);  % Both significant
                    
                    for mon = monkey_included_for_analysis
                        m = find(monkey_included_for_loading==mon);
                        plot(delta_psycho(methods_of_select{ms,1}&select_monkey{m}, comp_pair{cp}(1)+(option-1)*3),...
                            delta_psycho(methods_of_select{ms,1}&select_monkey{m}, comp_pair{cp}(2)+(option-1)*3),...
                            monkey_shape{m}, 'Color', 'k',  'MarkerSize', 8, 'LineWidth', 0.75);
                        
                        plot(delta_psycho(methods_of_select{ms,1}&select_monkey{m}&Sig1,comp_pair{cp}(1)+(option-1)*3),...
                            delta_psycho(methods_of_select{ms,1}&select_monkey{m}&Sig1, comp_pair{cp}(2)+(option-1)*3),...
                            monkey_shape{m},'Color', 'k','MarkerFaceColor',colorcode(comp_pair{cp}(1),:), 'MarkerSize',8, 'LineWidth',0.75);
                        plot(delta_psycho(methods_of_select{ms,1}&select_monkey{m}&Sig2,comp_pair{cp}(1)+(option-1)*3),...
                            delta_psycho(methods_of_select{ms,1}&select_monkey{m}&Sig2, comp_pair{cp}(2)+(option-1)*3),...
                            monkey_shape{m},'Color', 'k','MarkerFaceColor',colorcode(comp_pair{cp}(2),:), 'MarkerSize',8, 'LineWidth',0.75);
                        plot(delta_psycho(methods_of_select{ms,1}&select_monkey{m}&Sig_both,comp_pair{cp}(1)+(option-1)*3),...
                            delta_psycho(methods_of_select{ms,1}&select_monkey{m}&Sig_both, comp_pair{cp}(2)+(option-1)*3),...
                            monkey_shape{m},'Color', 'k','MarkerFaceColor',[1,0.7,0], 'MarkerSize',8, 'LineWidth',0.75);    % Orange dots indicate cases with both significance
                    end
                    
                    % Line of unity
                    max_scale = max(max(delta_psycho_select{ms}(:, comp_pair{cp}+(option-1)*3)));
                    min_scale = min(min(delta_psycho_select{ms}(:, comp_pair{cp}+(option-1)*3)));
                    axis([min_scale-0.2 max_scale*1.05 min_scale-0.2 max_scale*1.05]);
                    plot([min_scale-0.2 max_scale*1.05], [min_scale-0.2 max_scale*1.05], 'k--');
                    
                    plot([0 0],[min_scale-0.2 max_scale*1.05], 'k--');
                    plot([min_scale-0.2 max_scale*1.05], [0 0], 'k--');
                    
                    nticks = 5;
                    xticks(roundn(linspace(min_scale,max_scale,nticks),-1));
                    yticks(roundn(linspace(min_scale,max_scale,nticks),-1));
                    
                    % Two-paired Wilcoxon test
                    pp = signrank(delta_psycho_select{ms}(:, comp_pair{cp}(1)+(option-1)*3),...
                        delta_psycho_select{ms}(:, comp_pair{cp}(2)+(option-1)*3));
                    
                    % Spearman and Pearson correlation
                    [r_spear, p_spear] = corr(delta_psycho_select{ms}(:, comp_pair{cp}(1)+(option-1)*3),...
                        delta_psycho_select{ms}(:, comp_pair{cp}(2)+(option-1)*3),'type', 'Spearman');
                    [r_pear, p_pear] = corr(delta_psycho_select{ms}(:, comp_pair{cp}(1)+(option-1)*3),...
                        delta_psycho_select{ms}(:, comp_pair{cp}(2)+(option-1)*3), 'type', 'Pearson');
                    
                    % Linear Fitting
                    [para, S] = polyfit(delta_psycho_select{ms}(:, comp_pair{cp}(1)+(option-1)*3),...
                        delta_psycho_select{ms}(:, comp_pair{cp}(2)+(option-1)*3),  1);
                    xx = min_scale : 0.1: max_scale;
                    Y = polyval(para,xx);
                    hl = plot(xx,Y,'k-','LineWidth',2);
                    
                    % Add text
                    txt_test = {[ 'n = ', num2str(min(sum(~isnan(delta_psycho_select{ms}(:,comp_pair{cp}+(option-1)*3))))) ],...
                        ['p-value: ', num2str(roundn(pp,-3))]};
                    str_test = text(min_scale-0.2+0.05,max_scale, txt_test, 'HorizontalAlignment','left');
                    set(str_test, 'Color', [pp<0.05,(pp<0.05)*0.7,0]);
                    
                    txt_corr_spear = {['r_{Spearman} = ', num2str(roundn(r_spear,-3))], ['p_{Spearman}: ', num2str(roundn(p_spear,-3))]};
                    txt_corr_pear = {['r_{Pearson} = ', num2str(roundn(r_pear,-3))], ['p_{Pearson}: ', num2str(roundn(p_pear,-3))]};
                    str_corr_spear = text(max_scale-0.8, min_scale+0.5, txt_corr_spear,'HorizontalAlignment','left');
                    str_corr_pear = text(max_scale-0.8, min_scale+0.65, txt_corr_pear,'HorizontalAlignment','left');
                    
                    % X-axis and Y-axis label
                    xlabel(stim_type_list{comp_pair{cp}(1)},'Color',colorcode(comp_pair{cp}(1),:));
                    ylabel(stim_type_list{comp_pair{cp}(2)},'Color',colorcode(comp_pair{cp}(2),:));
                    
                    % Annotate median value
                    text(delta_psycho_median{ms}(:, comp_pair{cp}(1)+(option-1)*3),min_scale+0.1,...
                        '\downarrow', 'Color', colorcode(comp_pair{cp}(1),:), 'FontSize', 15, 'FontWeight', 'bold');
                    text(min_scale+0.1, delta_psycho_median{ms}(:, comp_pair{cp}(2)+(option-1)*3),...
                        '\leftarrow', 'Color', colorcode(comp_pair{cp}(2),:), 'FontSize', 15, 'FontWeight', 'bold');
                    
                    % Title
                    if cp == 1
                        if option ==BIAS
                            title('\DeltaBias','Color', 'k', 'FontSize',20,'FontWeight','bold');
                        else
                            title('\DeltaThreshold','Color','k', 'FontSize',20,'FontWeight','bold');
                        end
                    end
                    
                    axis square;
                end
            end
            SetFigure();
        end
        
    end

%%
    function f1p4(debug)       % Normalized Bias comparison across conditions
        if debug; dbstack; keyboard; end
        
        % Normalize bias shift by dividing Threshold of non-ustim condition
        % in corresponding session
        % Xuefei Yu, 2018, Neuron
        norm_bias = delta_psycho(:,1:3) ./ psycho_ctrl(:, 4:6);
                
        % Monkey mask has been included
        methods_of_select = {
            select_all , 'All sites'
            select_tsites , 'Typical sites'
            select_no_tsites, 'Non-typical sites'
            };
        
        comp_pair = {[1 2]; [1 3]; [2 3]};   % comparing pairs
        
        % Significance indication of Bias
%         sig_ind = p_value_psycho(:,1:3) <= 0.05;
        sig_ind = bootstrap_p_value(:,1:3) <= 0.01;
        
        
        %% 3D plotting for simultaneously comparison 
        % Added by ZZ @20231030
        set(figure, 'name','Normalized delta bias correlation'); hold on; 
        mon = find(monkey_included_for_loading == monkey_included_for_analysis);
        plot3(norm_bias(select_all,1), norm_bias(select_all,2), norm_bias(select_all,3),...
            monkey_shape{mon}, 'color', 'k','MarkerSize', 8, 'LineWidth', 0.75);
        
        % Significance indication 
        sig_condi_num = sum(sig_ind,2); 
        for c = 1:3
            plot3(norm_bias(select_all&sig_ind(:,c)&(sig_condi_num==1),1), norm_bias(select_all&sig_ind(:,c)&(sig_condi_num==1),2),...
                norm_bias(select_all&sig_ind(:,c)&(sig_condi_num==1),3),...
                monkey_shape{mon}, 'color','k','MarkerFaceColor',colorcode(c,:),'MarkerSize', 8, 'LineWidth', 0.75);
        end
        
        plot3(norm_bias(select_all&(sig_condi_num>1),1), norm_bias(select_all&(sig_condi_num>1),2), norm_bias(select_all&(sig_condi_num>1),3),...
            monkey_shape{mon}, 'color', 'k','MarkerFaceColor',[0.9290 0.6940 0.1250],'MarkerSize', 8, 'LineWidth', 0.75);

        %fitting using pca
        [coeff, score, ~,~,explained,mu] = pca(norm_bias(select_all,:));
        % Direction of the first PC
        dirVect = coeff(:,1); PC1variance = explained(1) / sum(explained); 
        
        % Plotting the fitted line
        xxx = [min(score(:,1))-.2 max(score(:,1))+.2]; 
        endpts = [mu+xxx(1)*dirVect'; mu+xxx(2)*dirVect'];
        plot3(endpts(:,1), endpts(:,2), endpts(:,3), '-', 'color', [0.9290 0.6940 0.1250], 'linew', 3);
        title(['The first PC explains ' num2str(PC1variance) ' variance']); 
        
                
        xlabel('Vestibular'); ylabel('Visual'); zlabel('Combined'); 
        grid on; 
        view(45, 25);
        SetFigure(); 

        
        for ms = 1:size(methods_of_select, 1)
            % Median value
            norm_delta_bias_median{ms} = nanmedian(norm_bias(methods_of_select{ms,1},:));
            % Selected value
            norm_delta_bias_select{ms} = norm_bias(methods_of_select{ms,1},:);
            
            % Kruskal-Wallis test (non-parametric one-way ANOVA) for stim types
            % With multiple comparions of Fisher's least significant difference procedure
            [p_stimtype{ms}, ~, stats] = kruskalwallis(norm_bias(methods_of_select{ms,1},:), stim_type_list, 'off');
            comp_result{ms} = multcompare(stats, 'CType','lsd', 'Display','off');
            
            % -- Plotting
            set(figure(figN), 'name', ['Across-condition comparison of Norm-Bias ' methods_of_select{ms, 2}], 'Position',[80 10 1000 1000]); clf; figN=figN+1;
            ha = tight_subplot(3,1,[.1 .05],[.1 .1], .1);
            
            for cp = 1:length(comp_pair)
                axes(ha(cp)); hold on;
                
                Sig1 = sig_ind(:, comp_pair{cp}(1)) & ~sig_ind(:, comp_pair{cp}(2));   % Only significant in the first condition
                Sig2 = sig_ind(:, comp_pair{cp}(2)) & ~sig_ind(:, comp_pair{cp}(1));   % Only significant in the second condition
                Sig_both = sig_ind(:, comp_pair{cp}(2)) & sig_ind(:, comp_pair{cp}(1));  % Both significant
                
                for mon = monkey_included_for_analysis
                    m = find(monkey_included_for_loading==mon);
                    plot(norm_bias(methods_of_select{ms,1}&select_monkey{m}, comp_pair{cp}(1)),...
                        norm_bias(methods_of_select{ms,1}&select_monkey{m}, comp_pair{cp}(2)),...
                        monkey_shape{m}, 'Color', 'k',  'MarkerSize', 8, 'LineWidth', 0.75);
                    
                    plot(norm_bias(methods_of_select{ms,1}&select_monkey{m}&Sig1,comp_pair{cp}(1)),...
                        norm_bias(methods_of_select{ms,1}&select_monkey{m}&Sig1, comp_pair{cp}(2)),...
                        monkey_shape{m},'Color', 'k','MarkerFaceColor',colorcode(comp_pair{cp}(1),:), 'MarkerSize',8, 'LineWidth',0.75);
                    plot(norm_bias(methods_of_select{ms,1}&select_monkey{m}&Sig2,comp_pair{cp}(1)),...
                        norm_bias(methods_of_select{ms,1}&select_monkey{m}&Sig2, comp_pair{cp}(2)),...
                        monkey_shape{m},'Color', 'k','MarkerFaceColor',colorcode(comp_pair{cp}(2),:), 'MarkerSize',8, 'LineWidth',0.75);
                    plot(norm_bias(methods_of_select{ms,1}&select_monkey{m}&Sig_both,comp_pair{cp}(1)),...
                        norm_bias(methods_of_select{ms,1}&select_monkey{m}&Sig_both, comp_pair{cp}(2)),...
                        monkey_shape{m},'Color', 'k','MarkerFaceColor',[1,0.7,0], 'MarkerSize',8, 'LineWidth',0.75);    % Orange dots indicate cases with both significance
                end
                
                % Line of unity
                max_scale = max(max(norm_delta_bias_select{ms}(:, comp_pair{cp})));
                min_scale = min(min(norm_delta_bias_select{ms}(:, comp_pair{cp})));
                axis([min_scale-0.2 max_scale*1.05 min_scale-0.2 max_scale*1.05]);
                plot([min_scale-0.2 max_scale*1.05], [min_scale-0.2 max_scale*1.05], 'k--');
                
                plot([0 0],[min_scale-0.2 max_scale*1.05], 'k--');
                plot([min_scale-0.2 max_scale*1.05], [0 0], 'k--');
                
%                 nticks = 5;
%                 xticks(roundn(linspace(min_scale,max_scale,nticks),-1));
%                 yticks(roundn(linspace(min_scale,max_scale,nticks),-1));
                
                % Change to multiple comparison kruskal-wallis test
                %                 % Two-paired Wilcoxon test
                %                 pp = signrank(norm_delta_bias_select{ms}(:, comp_pair{cp}(1)),...
                %                     norm_delta_bias_select{ms}(:, comp_pair{cp}(2)));
                
                % Spearman and Pearson correlation
                [r_spear, p_spear] = corr(norm_delta_bias_select{ms}(:, comp_pair{cp}(1)),...
                    norm_delta_bias_select{ms}(:, comp_pair{cp}(2)),'type', 'Spearman');
                [r_pear, p_pear] = corr(norm_delta_bias_select{ms}(:, comp_pair{cp}(1)),...
                    norm_delta_bias_select{ms}(:, comp_pair{cp}(2)), 'type', 'Pearson');
                
%                 % Linear Fitting
%                 [para, S] = polyfit(norm_delta_bias_select{ms}(:, comp_pair{cp}(1)),...
%                     norm_delta_bias_select{ms}(:, comp_pair{cp}(2)),  1);
                
                % Change to Perpendicular regression (PCA)
                fitType = regress_perp(norm_delta_bias_select{ms}(:, comp_pair{cp}(1)),...
                    norm_delta_bias_select{ms}(:, comp_pair{cp}(2)),1,0.05); 
                linPara(1) = fitType.k; linPara(2) = fitType.b; 
                
                xx = min_scale : 0.1: max_scale;
                Y = linPara(1)*xx + linPara(2); 
                
                % Constrain the line
                yy = norm_delta_bias_select{ms}(:, comp_pair{cp}(2));
                xx = xx(min(yy)<=Y & Y<=max(yy));
                Y = Y(min(yy)<=Y & Y<=max(yy)); 
%                 Y = polyval(para,xx);
                hl = plot(xx,Y,'k-','LineWidth',2);
                
                % Add text
                txt_test = {[ 'n = ', num2str(min(sum(~isnan(norm_delta_bias_select{ms}(:,comp_pair{cp}))))) ],...
                    ['p-value: ', num2str(roundn(comp_result{ms}(cp,end),-3))]};
                str_test = text(min_scale-0.2+0.05,max_scale, txt_test, 'HorizontalAlignment','left');
                set(str_test, 'Color', [comp_result{ms}(cp,end)<0.05,(comp_result{ms}(cp,end)<0.05)*0.7,0]);
                
                txt_corr_spear = {['r_{Spearman} = ', num2str(roundn(r_spear,-3))], ['p_{Spearman}: ', num2str(roundn(p_spear,-3))]};
                txt_corr_pear = {['r_{Pearson} = ', num2str(roundn(r_pear,-3))], ['p_{Pearson}: ', num2str(roundn(p_pear,-3))]};
                str_corr_spear = text(max_scale-0.8, min_scale+0.5, txt_corr_spear,'HorizontalAlignment','left');
                str_corr_pear = text(max_scale-0.8, min_scale+0.65, txt_corr_pear,'HorizontalAlignment','left');
                
                % X-axis and Y-axis label
                xlabel(stim_type_list{comp_pair{cp}(1)},'Color',colorcode(comp_pair{cp}(1),:));
                ylabel(stim_type_list{comp_pair{cp}(2)},'Color',colorcode(comp_pair{cp}(2),:));
                
                % Annotate median value
                plot(norm_delta_bias_median{ms}(:, comp_pair{cp}(1)),min_scale+0.1,...
                    'v', 'Color', colorcode(comp_pair{cp}(1),:), 'MarkerSize', 15);
                text(norm_delta_bias_median{ms}(:, comp_pair{cp}(1)),min_scale+0.1,...
                    num2str(norm_delta_bias_median{ms}(:, comp_pair{cp}(1))), 'Color', colorcode(comp_pair{cp}(1),:))
                plot(min_scale+0.1, norm_delta_bias_median{ms}(:, comp_pair{cp}(2)),...
                    '<', 'Color', colorcode(comp_pair{cp}(2),:), 'MarkerSize', 15);
                text(min_scale+0.1, norm_delta_bias_median{ms}(:, comp_pair{cp}(2)),...
                    num2str(norm_delta_bias_median{ms}(:, comp_pair{cp}(2))), 'Color', colorcode(comp_pair{cp}(2),:))
                
                % Title
                if cp == 1
                    %                     title('\DeltaBias','Color', 'k', 'FontSize',20,'FontWeight','bold');
                    title(['p_{kruskalwallis}: ' num2str(roundn(p_stimtype{ms},-3))]);
                end
                
                axis square;
            end
            SetFigure();
        end
    end

%%
    function f1p5(debug)       % Correct rate effct in each headings
        if debug; dbstack; keyboard; end
        
        % Monkey mask has been included
        methods_of_select = {
            select_all , 'All sites'
            select_tsites , 'Typical sites'
            select_no_tsites, 'Non-typical sites'
            };
        
        for ms = 1:size(methods_of_select,1)
            set(figure(figN), 'pos',[50 60 1000 500], 'name',['Delta Correct Rate - ' methods_of_select{ms,2}]); clf; hold on; figN=figN+1;
            
%             delta_CorrectRate{ms} = correct_rate_ustim(methods_of_select{ms,1},:,:) - correct_rate_ctrl(methods_of_select{ms,1},:,:);   % absolute value to control the different biased side
            delta_CorrectRate{ms} = abs(correct_rate_ustim(methods_of_select{ms,1},:,:) - correct_rate_ctrl(methods_of_select{ms,1},:,:));   % absolute value to control the different biased side
            delta_CorrectRate_mean{ms} = squeeze(mean(delta_CorrectRate{ms}));
            delta_CorrecRate_sem{ms} = squeeze(std(delta_CorrectRate{ms},0,1) / sqrt(size(delta_CorrectRate{ms},1)));
            
            % Delta correct rate in each heading
            for c = 1:3   % conditions
                % One-way Anova for headings
                p_heading(c) = anova1(squeeze(delta_CorrectRate{ms}(:,:,c)), unique_heading', 'off');
                
                for hh = 1:length(unique_heading)
                    % t-test
                    [~, p_delta_correctrate(hh,c)] = ttest(delta_CorrectRate{ms}(:,hh,c));
                end
            end
            
            % Plotting
            errorbar(repmat(unique_heading, 1, 3), delta_CorrectRate_mean{ms}, delta_CorrecRate_sem{ms}, 'LineWidth',2);
            xlim([-7 7]); ylims= ylim;
            plot([-9 9], [0 0], 'k--', 'linew',1);
            xlabel('Heading from ipsi to contra (o)'); ylabel('\Delta Correct Rate');
            legend(stim_type_list, 'AutoUpdate','off');
            
            % Significance indication
            sig_ind = p_delta_correctrate <=0.05;
            for c = 1:3
                if any(sig_ind(:,c))
                    text(unique_heading(sig_ind(:,c)), repmat(ylims(2)+(c-1)*0.005, sum(sig_ind(:,c)), 1),...
                        '\ast', 'color', colorcode(c,:), 'FontWeight','bold');
                end
                
                text(6, ylims(1)+c*0.01, ['p_{ANOVA}: ' num2str(p_heading(c))], 'color',colorcode(c,:));
            end
            
            title(['n = ' num2str(size(delta_CorrectRate{ms}, 1))]);
            
            SetFigure();
        end
    end

%%
    function f1p6(debug)       % Temporal dynamics of dPSE
        if debug; dbstack; keyboard; end
        
        % Monkey mask has been included
        methods_of_select = {
            select_all , 'All sites'
            select_tsites , 'Typical sites'
            select_no_tsites, 'Non-typical sites'
            };
        
        for ms = 1:size(methods_of_select,1)
            set(figure(figN), 'pos',[50 60 1000 500], 'name',['Delta Bias Dynamic - ' methods_of_select{ms,2}]); clf; hold on; figN=figN+1;
            
            dPSE_dynamic_mean{ms} = squeeze(nanmean(dPSE_dynamic(methods_of_select{ms,1},:,:)));
            dPSE_dynamic_sem{ms} = squeeze(nanstd(dPSE_dynamic(methods_of_select{ms,1},:,:), 0,1)) /...
                sqrt(sum(~isnan(dPSE_dynamic(methods_of_select{ms,1},1,1))));
            
            for c = 1:3
                % one-way ANOVA across different time window
                p_windows(c) = anova1(squeeze(dPSE_dynamic(methods_of_select{ms,1},:,c)),[], 'off');
                
                % t-test for each time window
                for nw = 1:shift_n   % number of window
                    [~, p_each_window(nw,c)] = ttest(dPSE_dynamic(methods_of_select{ms,1},nw,c));
                end
                
            % == Plotting
            xx = 1:shift_n;
            errorbar(xx', dPSE_dynamic_mean{ms}(:,c), dPSE_dynamic_sem{ms}(:,c),'o',...
                'MarkerSize',10, 'MarkerFaceColor',colorcode(c,:), 'linew',2);
                
            end
            
            
            xlim([0 shift_n+0.5]); ylims=ylim;
            plot(xlim, [0 0], 'k--', 'linew',1);
            
            sig_ind = p_each_window <= 0.01;
            for c = 1:3
                if any(sig_ind(:,c))
                    text(xx(sig_ind(:,c)), repmat(ylims(2)+(c-1)*0.01, 1,sum(sig_ind(:,c))),...
                        '\ast', 'color', colorcode(c,:), 'FontWeight','bold');
                end
                
                text(6, ylims(1)+c*0.01, ['p_{ANOVA}: ' num2str(p_windows(c))], 'color',colorcode(c,:));
            end
            
            xlabel('Temporal period'); ylabel('\Delta PSE (o)');
            SetFigure();
            
        end
        
    end


%% =============    For Correlations  ========================
% Sign DDI according to preferred direction relative to hemisphere
% Positive indicate contralateral preferrence 
signed_DDI = memsac_DDI .* ((memsac_DDI_if_contra - 0.5)*2); 

% the first 3 columns from uStim, delta PSE in three conditions
% MemSac_DDI, 6 columns
% HD_ChoiceDiv, 3 columns
ustim_memsac_HD_data = [dPSE signed_DDI HD_ChDiv];

%%
    function f2p1(debug)
        if debug; dbstack; keyboard; end
        
        % Monkey mask has been included
        methods_of_select = {
            select_all , 'All sites'
            };
        
        % Normalize the data to [-1, 1]
        nor_dPSE = mapminmax(dPSE(methods_of_select{1,1},:)')';
        nor_DDI = mapminmax(signed_DDI(methods_of_select{1,1},:)')';
        nor_ChDiv = mapminmax(HD_ChDiv(methods_of_select{1,1},:)')';
        ustim_MS_HD_sorted_normal = [nor_dPSE nor_DDI nor_ChDiv];
        
        % Because I doesn't run MemSac or Heading task ahead of Microstim
        % So there are some NaNs
        % For the beauty of the figure, I sort the sites according to ChoiceDive
        % firstly, and then if NaN in Heading, I sort by the order of MemSac_DDI
        nonNaN = ~isnan(nor_ChDiv(:,3));
        [~, ind1] = sort(nor_ChDiv(:,3)); % default, sorted by the choice preference in combined condition
        [~, ind2] = sort(nor_DDI(~nonNaN,3));   % secondary, sorted by the DDI during memory period
        
        nan_ind = ind1(end-length(ind2)+1:end, :);
        ind1(end-length(ind2)+1:end, :) = nan_ind(ind2);
        sort_ind(:,1) = ind1;
        
        % Second sorting method, just according to the value of dPSE
        [~, sort_ind(:,2)] = sort(nor_dPSE(:,1));    % ipsilateral (negaitve) first
        
        % == Heatmap Plotting
        % white = [255 255 255]/255; % for 0 updating
        % purple = [160 32 240]/255; % for positive updating
        % orange = [255 128 0]/255; % for negtive updating
        % color_data_neg = flipud((0:49)'/50 *(orange-white)+white);  color_data_posi = (0:49)'/50*(purple-white)+white;
        % color_map = [color_data_neg;color_data_posi];
        for i = 1:size(sort_ind,2)
            set(figure(figN),'Name','ustim-MemSac-HD_Heatmap'); clf; figN= figN+1;
            hm = heatmap(gcf,ustim_MS_HD_sorted_normal(sort_ind(:,i),:), 'colormap', spring,...
                'ColorLimits', [-1 1], 'ColorbarVisible', 'on', 'GridVisible', 'off', 'CellLabelColor','none', 'MissingDataColor','k');
            % hm.FontSize = 13;
            hm.YLabel = 'Cell Number';
            
            SetFigure();
        end
        
    end

%%
    function f2p2(debug)   % MemSac with dPSE Correlation scatter plotting
        if debug; dbstack; keyboard; end
        
        % Monkey mask has been included
        methods_of_select = {
            select_all , 'All sites'
            select_tsites , 'Typical sites'
            select_no_tsites, 'Non-typical sites'
            };
        
        %         ustim_memsac_HD_data = [dPSE signed_DDI HD_ChDiv];
        
        for ms = 1:size(methods_of_select,1)
            
%             % signed DDI with ustim results
%             h = LinearCorrelation(ustim_memsac_HD_data(methods_of_select{ms,1},6),...
%                 ustim_memsac_HD_data(methods_of_select{ms,1},1:3),...
%                 'EdgeColors',mat2cell(colorcode,[1;1;1]), 'FaceColors',{'none', 'none', 'none'}, 'Markers', {'o'},'MarkerSize', 10,...
%                 'LineStyles',{'b-','r-','g-'}, 'figN', figN, 'Xlabel', 'Aligned DDI(mem)', 'Ylabel', 'delta PSE'); figN = figN+1; hold on; 
%             
%             for c = 1:3
%                 plot(ustim_memsac_HD_data(methods_of_select{ms,1}&pPSE(:,c)<0.05,6),...
%                     ustim_memsac_HD_data(methods_of_select{ms,1}&pPSE(:,c)<0.05,c),...
%                     'o','color',colorcode(c,:),'MarkerSize', 10, 'MarkerFaceColor',colorcode(c,:));
%             end
%             SetFigure();
%            
            
            % Raw DDI with abs(dPSE)
            h = LinearCorrelation(memsac_DDI(methods_of_select{ms,1},3),...
                abs(dPSE(methods_of_select{ms,1},:)),...
                'EdgeColors',mat2cell(colorcode,[1;1;1]),'FaceColors',{'none', 'none', 'none'}, 'Markers', {'o'},'MarkerSize', 10,...
                'LineStyles',{'b-','r-','g-'}, 'figN', figN, 'Xlabel', 'DDI(mem)', 'Ylabel', 'abs(delta PSE)');  figN = figN+1;  hold on; 
            
            % Mark significant dPSE
            for c = 1:3
                %         plot(abs(ustim_memsac_HD_data(pPSE(:,1)<0.05,6)),abs(ustim_memsac_HD_data(pPSE(:,1)<0.05,1)), 'b.','MarkerSize', 30,'linew',2);
                %         plot(abs(ustim_memsac_HD_data(pPSE(:,2)<0.05,6)),abs(ustim_memsac_HD_data(pPSE(:,2)<0.05,2)), 'r.','MarkerSize', 30,'linew',2);
                plot(memsac_DDI(methods_of_select{ms,1}&bootstrap_p_bias(:,c)<0.01,6),...
                    abs(dPSE(methods_of_select{ms,1}&bootstrap_p_bias(:,c)<0.01,c)),...
                    'o','color',colorcode(c,:),'MarkerSize', 10,'MarkerFaceColor',colorcode(c,:));
            end
            
            SetFigure();
            
        end
    end

%%         
    function f2p3(debug)  % Heading ChoicePref Correlation with dPSE
        if debug; dbstack; keyboard; end
        
                % Monkey mask has been included
        methods_of_select = {
            select_all , 'All sites'
            select_tsites , 'Typical sites'
            select_no_tsites, 'Non-typical sites'
            };

        m = find(monkey_included_for_loading==monkey_included_for_analysis);
        
        for ms = 1: size(methods_of_select, 1)
            h = LinearCorrelation(HD_ChDiv(methods_of_select{ms,1},:),...
                dPSE(methods_of_select{ms,1},:),...
                'EdgeColors',mat2cell(colorcode,[1;1;1]),'FaceColors',{'none', 'none', 'none'}, 'Markers', monkey_shape(m),'MarkerSize', 10,...
                 'LineStyles',{'b-','r-','g-'}, 'figN', figN, 'Xlabel', 'ChoicePref', 'Ylabel', '\Delta PSE (o)');   figN = figN+1;  hold on; 
             
            % Mark significant dPSE
             for c = 1:3
                plot(HD_ChDiv(methods_of_select{ms,1}&bootstrap_p_bias(:,c)<0.01,c),...
                    dPSE(methods_of_select{ms,1}&bootstrap_p_bias(:,c)<0.01,c),...
                    monkey_shape{m},'color',colorcode(c,:),'MarkerSize', 10,'MarkerFaceColor',colorcode(c,:));
             end
            
        end
                
    end

%%  
    function f2p4(debug)  % Memsac DDI with Heading ChoicePref
        if debug; dbstack; keyboard; end
        
        % Monkey mask has been included
        methods_of_select = {
            select_all , 'All sites'
            select_tsites , 'Typical sites'
            select_no_tsites, 'Non-typical sites'
            };

        for ms = 1: size(methods_of_select,1)
            % DDI(mem) with abs(CDiv)
            h = LinearCorrelation(memsac_DDI(methods_of_select{ms,1},3),abs(HD_ChDiv(methods_of_select{ms,1},:)),...
                'XLabel','DDI(mem)', 'YLabel','abs(CDiv)',...
                'EdgeColors',mat2cell(colorcode,[1;1;1]),'FaceColors',{'none','none','none'},'Markers',{'o','o','o'},'MarkerSize', 10,...
                'LineStyles',{'b-','r-','g-'},'figN',figN);   figN = figN+1; hold on;
            
            for c = 1:3
                plot(memsac_DDI(methods_of_select{ms,1}&bootstrap_p_bias(:,c)<0.01,3),...
                    abs(HD_ChDiv(methods_of_select{ms,1}&bootstrap_p_bias(:,c)<0.01,c)),...
                    'o','color',colorcode(c,:),'MarkerSize', 10,'MarkerFaceColor',colorcode(c,:));
            end
            
            % Show individual cell selected from the figure
            h_line = plot(memsac_DDI(methods_of_select{ms,1},3), abs(HD_ChDiv(methods_of_select{ms,1},c)), 'visible','off'); hold on; 
            set([gca [h.group.dots] h_line], 'ButtonDownFcn',{@Show_individual_cell,h_line,select_all});
            
            SetFigure();
        end
                
    end


site_position = [];
cell_position();
%%
    function f3p1(debug)     % Microstimulation site distribution
        if debug; dbstack; keyboard; end
        
        % For specific monkey and hemisphere
        to_plot = {monkey_included_for_analysis; 1};
        
        select_mon_hems_ind = (site_position(:,1)==to_plot{1} & site_position(:,2)==to_plot{2});
        
        if to_plot{1} == 15 
            txt_title = 'Monkey D    Left Hemisphere';
        elseif to_plot{1} == 13
            txt_title = 'Monkey F    Left Hemisphere';
        end
        
        % Along the AP axis
        nbins = 10; 
        min_AP = min(site_position(select_mon_hems_ind,3));
        max_AP = max(site_position(select_mon_hems_ind,3));
        xedges_AP = linspace(floor(min_AP),ceil(max_AP),nbins);
        xbins_AP = (xedges_AP(1:end-1)+xedges_AP(2:end))/2;
        
        % Sites with significant dPSE in any stimulus conditions
        sig_site_ind = select_mon_hems_ind &  (any(bootstrap_p_bias<=0.01,2));
        sig_site_AP_counts = hist(site_position(sig_site_ind,3),xbins_AP);
        
        non_sig_site_ind = select_mon_hems_ind &(all(bootstrap_p_bias>0.01,2));
        non_sig_site_AP_counts = hist(site_position(non_sig_site_ind,3),xbins_AP);
        
        % Mann-Whitney U-test
        AP_sig = ranksum(site_position(sig_site_ind,3),site_position(non_sig_site_ind,3));
        depth_sig = ranksum(site_position(sig_site_ind,4),site_position(non_sig_site_ind,4));
        
        % AP Axis Plotting
        set(figure(figN),'pos',[80,50,1600,800], 'Name', txt_title); clf;   figN=figN+1;
        hs = subplot(1,2,1); hold on;
        h_bar = bar(xbins_AP,sig_site_AP_counts, 'histc');
        set(h_bar, 'EdgeColor','k',  'FaceColor','k');
        
        h_bar = bar(xbins_AP,-non_sig_site_AP_counts, 'histc');
        set(h_bar, 'EdgeColor','k',  'FaceColor',[0.5 0.5 0.5]);
        
        set(gca, 'xDir', 'reverse');    % Left -- Anterior
        title('Along AP Axis')
        xlabel({'AC','( mm )'});
        %         if AP_sig <= 0.05
        text(max(xbins_AP),max(sig_site_AP_counts),num2str(roundn(AP_sig,-4)));
        %         end
        
        % Depth Plotting
        min_depth = min(site_position(select_mon_hems_ind,4));
        max_depth = max(site_position(select_mon_hems_ind,4));
        xedges_depth = (floor(min_depth/500):ceil(max_depth/500))*500;
        xbins_depth = ((xedges_depth(1:end-1))+(xedges_depth(2:end)))/2;
        
        sig_site_depth_counts = hist(site_position(sig_site_ind,4), xbins_depth);
        non_sig_site_depth_counts = hist(site_position(non_sig_site_ind,4), xbins_depth);
        
        hs = subplot(1,2,2); hold on;
        h_bar = bar(xbins_depth,sig_site_depth_counts, 'histc');
        set(h_bar, 'EdgeColor','k',  'FaceColor','k');
        h_bar = bar(xbins_depth,-non_sig_site_depth_counts, 'histc');
        set(h_bar, 'EdgeColor','k',  'FaceColor',[0.5 0.5 0.5]);
        
        legend({'Any-Cell'; 'Non-Cell'}, 'Position',[0.9 0.9 0.1 0.1],'AutoUpdate','off');
        title('Along DV Axis');
        xlabel({'Depth into CD', '( \mum)'});
        %         if depth_sig <= 0.05
        text(0,max(sig_site_depth_counts),num2str(roundn(depth_sig,-4)));
        %         end
        
        set(gca, 'view', [-90,-90]);
        
        SetFigure();
    end

%%
    function f3p2(debug)     % Microstimulation site distribution comparison across monkeys
        if debug; dbstack; keyboard; end
        
        % For specific monkey and hemisphere
        to_plot = {monkey_included_for_analysis; 1};
        
        select_mon_hems_ind = (site_position(:,1)==to_plot{1} & site_position(:,2)==to_plot{2});
        sig_site_ind = select_mon_hems_ind &  (any(bootstrap_p_bias<=0.01,2));
        
        AP_loca_all = {site_position(select_mon_hems_ind(:,1),3); site_position(select_mon_hems_ind(:,2),3)}; 
        AP_loca_sig = {site_position(sig_site_ind(:,1),3); site_position(sig_site_ind(:,2),3)};
        
        % Mann-Whitney U-test
        % whether all stim location significant ?
        sig_all = ranksum(AP_loca_all{1}, AP_loca_all{2});
        % whether sig stim location significant ?
        sig_sig = ranksum(AP_loca_sig{1}, AP_loca_sig{2}); 
        
        % == Plotting
        nbins = 10; 
        % Along the AP axis
        min_AP = min(site_position(:,3));
        max_AP = max(site_position(:,3));
        xedges_AP = linspace(floor(min_AP),ceil(max_AP),nbins);
        xbins_AP = (xedges_AP(1:end-1)+xedges_AP(2:end))/2;
        
        % (1) All sites
        M1_all = hist(AP_loca_all{1}, xbins_AP);  % monkey 1
        M2_all = hist(AP_loca_all{2}, xbins_AP);  % monkey 2
        % (2) Sig sites
        M1_sig = hist(AP_loca_sig{1}, xbins_AP);  % monkey 1
        M2_sig = hist(AP_loca_sig{2}, xbins_AP);  % monkey 2
        
        
        set(figure(figN),'pos',[80,50,1600,800], 'Name', 'Site-distribution comparison'); clf;   figN=figN+1;
        
        % All sites
        hs = subplot(1,2,1); hold on;
        h_bar = bar(xbins_AP,M1_all, 'histc');
        set(h_bar, 'EdgeColor','k',  'FaceColor','k');
        
        h_bar = bar(xbins_AP,-M2_all, 'histc');
        set(h_bar, 'EdgeColor','k',  'FaceColor',[0.5 0.5 0.5]);
        
        set(gca, 'xDir', 'reverse');    % Left -- Anterior
        title('All Sites')
        xlabel({'AC','( mm )'});
        %         if AP_sig <= 0.05
        text(max(xbins_AP),max(M1_all),num2str(sig_all));
        %         end

        
        % Significant sites
        hs = subplot(1,2,2); hold on;
        h_bar = bar(xbins_AP,M1_sig, 'histc');
        set(h_bar, 'EdgeColor','k',  'FaceColor','k');
        h_bar = bar(xbins_AP,-M2_sig, 'histc');
        set(h_bar, 'EdgeColor','k',  'FaceColor',[0.5 0.5 0.5]);
        
                set(gca, 'xDir', 'reverse');    % Left -- Anterior

        legend({'M1'; 'M2'}, 'Position',[0.9 0.9 0.1 0.1],'AutoUpdate','off');
        title('Significant sites');
        %         if depth_sig <= 0.05
        text(0,max(M1_sig),num2str(sig_sig));
        %         end
        
%         set(gca, 'view', [-90,-90]);
        
        SetFigure();

    end
%% 
% Show individual cell selected from the figure
    function Show_individual_cell(~,~,h_line,select_for_this,couples)
        
        if nargin < 5
            couples = 1;  % show coupled cells in the figure
        end
        
        h_marker = guidata(gcbo);
        
        allX = get(h_line,'xData');
        allY = get(h_line,'yData');
        n_single_group = length(allX) / couples;
        
        % -------- Recover the cell number ------------
        if ismember('control', get(gcf,'currentModifier'))
            fileN = input('Which cell do you want from the figure?   ', 's'  );
            available_cells = find(select_for_this);
            
            if fileN(1) == '#'   % Direct input the orniginal cell #
                ori_cell_no = str2double(fileN(2:end));
                ind = sum(select_for_this(1:ori_cell_no));
                
            else
                nn = 1; iffind = [];
                while nn <= length(available_cells) && isempty(iffind)   % find the first match 
                    iffind = strfind(group_Microstim_result(available_cells(nn)).cellID{1}{1}, fileN);
                    nn = nn +1;
                end
                
                if ~isempty(iffind)
                    ind = nn - 1;
                else
                    fprintf('Are you kidding...?\n');
                    return
                end
                
                ori_cell_no = find(cumsum(select_for_this)==ind, 1); 
            end
            
        else   % Select cell from figure
            pos = get(gca, 'currentPoint'); posX = pos(1,1); posY = pos(1,2);
            [min_dis,ind] = min(abs(((posX-allX)/range(xlim)).^2+((posY-allY)/range(ylim)).^2));
            if min_dis > (range(xlim)^2+range(ylim)^2)/100 +inf ; return; end
            
            ind = mod(ind-1, n_single_group) + 1;   % Deal with coupled cells
            ori_cell_no = find(cumsum(select_for_this)==ind, 1);
        end
        
        % Plotting
        if ~isempty(h_marker); try delete(h_marker); catch; end; end
        
        all_inds = mod(1:length(allX), n_single_group) == mod(ind, n_single_group);  %Deal with coupled cells
        h_marker = plot(allX(all_inds), allY(all_inds), 'x', 'color','m', 'markersize',15, 'linew',3);
        
        % Do plot
        Plot_Microstim([], [], ori_cell_no);
        
        guidata(gcbo,h_marker);
        
    end

    function Plot_Microstim(~,~,ori_cell_no, h_1463_axes)
        
        
    end


end

%%  LIP HD Pooled Data
% Begin @ HH 201406
% Addaped by ZZ 2021

function function_handles = Group_HD(XlsData)
% Matlab version
matlab_ver = version;
ver_num = version('-release');
ver_num = str2num(ver_num(1:end-1));

%% Constants
LEFT = 1;
RIGHT = 2;
tmp1 = []; tmp2 = []; tmp3 = []; tmp4 = []; tmp5  = [];

%% Get data
num = XlsData.num;
txt = XlsData.txt;
raw = XlsData.raw;
header = XlsData.header;
stim_type_num = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======= Change the file directory===========

% Batch address
mat_address = {
    % Major protocol goes here (Address, Suffix)
    
    'D:\Paper_rawdata\Raw_data\CN\Recordings\CN_m15&13_SU_Heading', 'PSTH'
    
    % Associative protocols
    'D:\Paper_rawdata\Raw_data\CN\Recordings\CN_m15&13_SU_MemSac' , 'MemSac'
    
    };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask_all = {
    (strcmp(txt(:,header.Protocol),'HD') & (strcmpi(txt(:,header.Area),'CD'))) ...
    & (num(:,header.HD_rep) >= 8) & (num(:,header.Units_RealSU) == 1) & (num(:,header.Monkey) == 15 | num(:,header.Monkey) == 13 )...
    & (num(:,header.Chan1) < 10);  % To speed up loading. To include HD_dt. HH20160918
    
    (strcmp(txt(:,header.Protocol),'MemSac')) & (strcmpi(txt(:,header.Area),'CD')) ...
    & (num(:,header.Chan1) < 10) & (num(:,header.Monkey) == 15 | num(:,header.Monkey) == 13);
    }; % Now no constraint on monkeys


% Add flexible monkey mask here (but I've still decided to choose monkey for analysis below). HH20150723
monkey_included_for_loading = [15 13];
monkey_included_for_analysis = [15 13];
monkey_marker = {'o', '^'};
monkey_line = {'-', '--'};

monkey_mask_for_loading = false(size(num,1),1);
for mm = 1:length(monkey_included_for_loading)
    monkey_mask_for_loading = monkey_mask_for_loading | (num(:,header.Monkey) == monkey_included_for_loading(mm));
end

% Now apply monkey mask
for mm = 1:length(mask_all)
    mask_all{mm} = mask_all{mm} & monkey_mask_for_loading;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for pp = 1:size(mat_address,1)
    % Xls Data
    xls_num{pp} = num(mask_all{pp},:);
    xls_txt{pp} = txt(mask_all{pp},:);
    xls_raw{pp} = raw(mask_all{pp},:);
    
    % Basic information : E.g.  '20140721m05s034h2x04y12d08478u6m5c67r2'
    cell_info_tmp = xls_num{pp}(:,[header.Date:header.Yloc header.Depth header.Chan1]);
    
    % Override MU with SU. HH20150422
    cell_info_tmp(xls_num{pp}(:,header.Units_RealSU)==1 & xls_num{pp}(:,header.Chan1)==1,end) = 5;
    
    cell_info{pp} = strsplit(sprintf('%8dm%02ds%03dh%01dx%02dy%02dd%05du%02d\n',cell_info_tmp'),'\n')';
    cell_info{pp} = strcat(cell_info{pp}(1:end-1),xls_txt{pp}(:,header.FileNo));
end
cd(mat_address{1,1});

%% Establish mat-to-xls relationship and load data into group_result.mat
% %{

global group_result; % This lets me reload the whole m.file without loading group_result, which speeds up my debugging.

%

if isempty(group_result)
    group_result(size(xls_txt{1},1)).cellID = [];  % The unique cellID
    load_group_result();
end

    function load_group_result()
        %
        % Load .mat files and put the data into a large structure array "group_result(i).mat_raw"
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tolerance_enable = 1;
        depthTol = 0; % Depth tolerance
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tic;
        progressbar('Load .mat files');
        not_match = 0;
        
        for major_i = 1:length(group_result) % Major protocol loop
            for pp = 1:size(mat_address,1)
                if pp == 1 % Major protocol
                    % Get basic information for cellID
                    file_i = major_i;
                else % Search for corresponding files in xls
                    match_i = find(strncmp(cell_info{1}{major_i},cell_info{pp},32));
                    
                    if length(match_i) == 1 % Exact match
                        file_i = match_i;
                    elseif length(match_i) > 1
                        
                        % Deal with 2point memsac. @HH20160906
                        if_2pt_memsac = strcmp('2MemSac',xls_txt{pp}(match_i,header.Protocol));
                        
                        if sum(if_2pt_memsac)  % We have 2pt memsac for this cell (almost all Messi cells but none for Polo)
                            % Note now the file_i has a length of 2.
                            % fprintf('Memsac and 2p-Memsac detected\n');
                            file_i(1) = match_i(find(if_2pt_memsac==0,1,'first')); % The first 8pt Memsac
                            file_i(2) = match_i(find(if_2pt_memsac==1,1,'first')); % The first 2pt Memsac
                        else  % No 2pt memsac, real duplication (maybe different eccentricities)
                            fprintf('More than one match have been found for %s,\n Major ID = %s\n', xls_txt{1}{major_i,header.FileNo},cell_info{1}{major_i});
                            disp(cell_info{pp}(match_i));
                            %                 file_choose = input('Which one do you want?');
                            %                 file_i = match_i(file_choose);
                            file_i = match_i(1);  % The first appearance by default.
                        end
                        
                        %                 keyboard;
                    else % Special cases for inexact match (MU-SU)
                        if ~tolerance_enable || ~strcmp(xls_txt{1}(major_i,header.Note),'MU=SU')
                            fprintf('No exact matched files for %s, ID = %s\n', xls_txt{1}{major_i,header.FileNo},cell_info{1}{major_i});
                            not_match = not_match + 1;
                            file_i = NaN;
                        else  % Unit tolerance enabled and "MU=SU"   HH20160920
                            match_i = find(strncmp(cell_info{1}{major_i},cell_info{pp},29));
                            
                            if length(match_i) == 1
                                fprintf('Unit tolerance found:\n   %s (major)\n   %s\n\n', cell_info{1}{major_i}, cell_info{pp}{match_i});
                                file_i = match_i;
                            elseif length(match_i) > 1
                                % Deal with 2point memsac. @HH20160906
                                if_2pt_memsac = strcmp('2MemSac',xls_txt{pp}(match_i,header.Protocol));
                                
                                if sum(if_2pt_memsac)  % We have 2pt memsac for this cell (almost all Messi cells but none for Polo)
                                    % Note now the file_i has a length of 2.
                                    % fprintf('Memsac and 2p-Memsac detected\n');
                                    file_i(1) = match_i(find(if_2pt_memsac==0,1,'first')); % The first 8pt Memsac
                                    file_i(2) = match_i(find(if_2pt_memsac==1,1,'first')); % The first 2pt Memsac
                                    fprintf('Unit tolerance found (and 2pt memsac):\n   %s (major)\n   %s\n\n   %s\n\n', cell_info{1}{major_i}, cell_info{pp}{match_i(1)},cell_info{pp}{match_i(2)});
                                else  % No 2pt memsac, real duplication (maybe different eccentricities)
                                    fprintf('More than one match have been found for %s,\n Major ID = %s\n', xls_txt{1}{major_i,header.FileNo},cell_info{1}{major_i});
                                    disp(cell_info{pp}(match_i));
                                    %                 file_choose = input('Which one do you want?');
                                    %                 file_i = match_i(file_choose);
                                    file_i = match_i(1);  % The first appearance by default.
                                end
                                
                            else
                                % Depth tolerance
                                match_i = find(strncmp(cell_info{1}{major_i},cell_info{pp},23),1);
                                depDiff = abs(xls_num{pp}(match_i,header.Depth)-xls_num{1}(major_i,header.Depth));
                                
                                if any(depDiff <= depthTol) % Within depth tolerance
                                    match_i = match_i(find(depDiff == min(depDiff),1));
                                    fprintf('Depth tolerance found:\n   %s (major)\n   %s\n\n', cell_info{1}{major_i}, cell_info{pp}{match_i});
                                    file_i = match_i;
                                else
                                    fprintf('!! All tolerance failed:\n    %s (major)\n\n',cell_info{1}{major_i});
                                    file_i = NaN;
                                    not_match = not_match + 1;
                                end
                            end
                        end
                        
                    end
                end
                
                % Load .mat for major and associative protocols
                if ~isnan(file_i)
                    try
                        for ii = 1:length(file_i) % For 2pt_memsac. @HH20160906
                            mat_file_name = sprintf('%s_%g',xls_txt{pp}{file_i(ii),header.FileNo},xls_num{pp}(file_i(ii),header.Chan1));
                            mat_file_fullname = [mat_address{pp,1} '\' mat_file_name '_' mat_address{pp,2}];
                            
                            raw = load(mat_file_fullname);
                            
                            group_result(major_i).cellID{pp}{ii} = cell_info{pp}{file_i(ii)};
                            
                            if strcmp('2MemSac',xls_txt{pp}(file_i(ii),header.Protocol)) % 2Mem-sac
                                group_result(major_i).(['mat_raw_2pt' mat_address{pp,2}]) = raw.result;  % Dynamic structure
                            else
                                group_result(major_i).(['mat_raw_' mat_address{pp,2}]) = raw.result;  % Dynamic structure
                            end
                        end
                    catch err
                        fprintf('Error Loading %s\n',[mat_file_name '_' mat_address{pp,2}]);
                        keyboard;
                    end
                else
                    group_result(major_i).(['mat_raw_' mat_address{pp,2}]) = [];
                end
                
            end
            
            progressbar(major_i/length(group_result));
            
        end
        fprintf('Loaded %g files (%g of them are incomplete).\n',length(group_result),not_match);
        toc
        
    end

%% Get some values from group_result(i).mat_raw_xxx our to group_result(i) for easier access
%  Preprocess each .mat_raw for different purposes (which have not been done in Batch Processing)

% HH20150414: Modality divergence added.
% HH20150419: Choice Divergence/Preference and Modality Divergence/Preference have been moved to Batch processing
% HH20150419: Different temporal alignment have been included

%%%%%%% Order corresponding to "Sort_Id" in TEMPO_GUI processing %%%%
ALL_CorrectCHOICE = 1; CORRECT_ANGLE = 2; CHOICE_DIFFICULT = 3; OUTCOME = 4; WRONG_ANGLE = 5; CORRECTNWRONG_ANGLE = 6; All_Choice = 7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
representative_cell = 94; % Last representative cell (HD, not HD_dt) for Duncan. ZZ20201207
% counts = 0;
for i = 1:length(group_result)
    
    group_result(i).repN = group_result(i).mat_raw_PSTH.repetitionN;
    unique_stim_type = group_result(i).mat_raw_PSTH.unique_stim_type;  % In case we don't have all three conditions
    group_result(i).length_unique_stim_type = length(unique_stim_type);
    group_result(i).unique_heading = group_result(i).mat_raw_PSTH.CP{1,1}.raw_CP_result{1}.Neu_tuning(:,1);
    
    % ------------ Patch in HD_dt cells. HH20160919 -----------
    % To include 18 HD_dt typical cells into my original HD tasks, I align PSTH using the center of the Gaussian movement instead of the
    % actual stim-on for those cells. So the time_alignment{1} should be moved 125 ms forward (in HD_dt task, stim duration = 1750 ms)
    % So now, VSTIM_ON_CD   means --> 1500 ms Gaussian begins (125 ms after the real VSTIM_ON_CD)
    %         VSTIM_OFF_CD  means --> 1500 ms Gaussian ends (125 ms before the real VSTIM_OFF_CD)
    % Here I just mark these cells out.
    if isfield(group_result(i).mat_raw_PSTH,'HD_dt_patched')
        group_result(i).HD_dt_patched = 1;
    else
        group_result(i).HD_dt_patched = 0;
    end
    
    group_result(i).time_marker{1} = [0 mean(group_result(i).mat_raw_PSTH.align_offsets_others{1})];
    group_result(i).time_marker{2} = [0];
    group_result(i).time_marker{3} = [0];
    
    % Psychophysics
    group_result(i).Psy_para = nan(2,3);
    for stim_type = 1:3
        k = find(stim_type == unique_stim_type);
        if ~isempty(k)   % We have this condition
            group_result(i).Psy_para(1,stim_type) = group_result(i).mat_raw_PSTH.CP{1,k}.Psy_para(2); % Threshold
            group_result(i).Psy_para(2,stim_type) = group_result(i).mat_raw_PSTH.CP{1,k}.Psy_para(1); % Bias
        end
    end
    
    %     for j = 1:2 % For two temporal alignments
    for j = 1:3 % For Three temporal alignments  % changed by ZZ 20210506
        
        
        % 1). For CP, here I decide to use pref. direction defined in Anne 2014 paper (100-200
        % ms before decision). But in the previous CP_HH.m, I assigned pref. direction for each time bin.
        % So here I check whether the local pref. direction is aligned with the "mode" (Zhong Shu) preferred
        % direction several time bins before decision. If not, I flip the original CP value.
        %  XXX @HH20141203 NO NEED To FLIP ANYMORE (I did this in HD_PSTH_HH)
        
        CP_ts = group_result(i).mat_raw_PSTH.CP{j,1}.ts;
        group_result(i).CP_ts{j} = CP_ts;
        
        
        % *** NOTE: I have done "always output three conditions" in TEMPO_GUI for all data except CP and PSTH. @HH20150419
        group_result(i).CP{j} = nan(3,length(CP_ts));
        group_result(i).CP_p_value{j} = nan(3,length(CP_ts));
        
        for stim_type = 1:3  % Always output three conditions
            
            k = find(stim_type == unique_stim_type);
            
            if ~isempty(k)   % We have this condition
                
                group_result(i).CP{j}(stim_type,:) = group_result(i).mat_raw_PSTH.CP{j,k}.CP_grand;
                
                % CP p value (permutation 1000)
                group_result(i).CP_p_value{j}(stim_type,:) = group_result(i).mat_raw_PSTH.CP{j,k}.CP_p;
            end
        end
        
        % 2). Rate p_value is straightforward
        
        rate_ts_temp = group_result(i).mat_raw_PSTH.PSTH{j,ALL_CorrectCHOICE,1}.ts;
        group_result(i).rate_ts{j} =  rate_ts_temp;
        
        group_result(i).rate_p_value{j} = NaN(3,length(rate_ts_temp));
        
        for stim_type = 1:3  % Always output three conditions
            
            k = find(stim_type == unique_stim_type);
            if ~isempty(k)   % We have this condition
                group_result(i).rate_p_value{j}(stim_type,:) =  group_result(i).mat_raw_PSTH.PSTH{j,ALL_CorrectCHOICE,1}.ps(k,:);
            end
            
        end
    end
    
    group_result(i).PREF_PSTH = group_result(i).mat_raw_PSTH.PREF;
    
    % 3). Next I calculate the "Choice divergence": AUC between different
    % choices under different conditions and difficulty levels.
    
    group_result(i).ChoicePreference = group_result(i).mat_raw_PSTH.ChoicePreference;
    
    % In TEMPO_GUI, choice preference uses the cell's PREF as its preferred direction
    % (because the hemisphere is unknown unless accessible to Result.xls)
    % Now I transform it to be related to "Contralateral" (Anne 2014)
    group_result(i).if_contralateral = xls_num{1}(i,header.Hemisphere) ~= group_result(i).PREF_PSTH;
    group_result(i).ChoicePreference = group_result(i).ChoicePreference * sign(group_result(i).if_contralateral - 0.5);
    
    group_result(i).ChoicePreference_pvalue =  group_result(i).mat_raw_PSTH.ChoicePreference_pvalue;
        
    % 4). Here comes the "modality divergence" (1-2,1-3,2-3)
    %     AUC between different modalities (regardless of choice)     HH20150415
    
    group_result(i).ModalityPreference = group_result(i).mat_raw_PSTH.ModalityPreference;
    group_result(i).ModalityPreference_pvalue =  group_result(i).mat_raw_PSTH.ModalityPreference_pvalue;
    
    
    % 5) Mem-sac dynamics (from .mat file) goes here
    if ~isempty(group_result(i).mat_raw_MemSac)
        group_result(i).MemSac_p = group_result(i).mat_raw_MemSac.p;
        group_result(i).MemSac_DDI = group_result(i).mat_raw_MemSac.DDI;
        group_result(i).MemSac_vectSum = group_result(i).mat_raw_MemSac.vectSum;
        group_result(i).MemSac_AI = group_result(i).mat_raw_MemSac.activityIndex;
        
        % Only align to saccade here
        group_result(i).MemSac_ts = group_result(i).mat_raw_MemSac.t_centers{3};
        group_result(i).MemSac_PSTH = reshape([group_result(i).mat_raw_MemSac.result_PSTH_anne_mean{3,:}],length(group_result(i).MemSac_ts),[]);
        
    end
    
    
    % 6) Spike waveform (finally!). HH20160920
    waveform = str2num(cell2mat(xls_txt{1}(i,header.meanWav)));
    % waveform = smooth(waveform,5)'; waveform = waveform - min(waveform); waveform = waveform./max(waveform);
    % Because I recorded using CED (sampling rate: 25kHz) at the beginning,
    % and then AO (22kHz)
    % ZZ @20231011
    ci = strfind(group_result(i).mat_raw_PSTH.FILE, 'c');   % cell ind
    ri = strfind(group_result(i).mat_raw_PSTH.FILE, 'r');   % run ind
    cell_num = str2num(group_result(i).mat_raw_PSTH.FILE(ci+1: ri-1));  
    group_result(i).cell_num = cell_num; 
    if cell_num < 200
         dt = 1/25; % 1k ms / 25 kHz
    else
        dt = 1/22;
    end
    interp_dt = 1/1000; % Spline interp to 1us. (Keven Johnston 2009 JNS)
            
    if ~isempty(waveform)
        
        waveform = spline(0:dt:dt*(length(waveform)-1),waveform,0:interp_dt:dt*(length(waveform)-1));
        
            waveform_t1 = find(waveform == max(waveform)); % Spike peak
        
        waveform_t_left = find(waveform == min(waveform(1:waveform_t1)),1); % Left trough
        waveform_t_right = find(waveform == min(waveform(waveform_t1:end)),1); % Right trough
        
        group_result(i).Waveform_width = waveform_t_right - waveform_t_left; 
        
        if cell_num < 200
            group_result(i).Waveform_peakToTrough = (waveform_t_right - waveform_t1)*interp_dt;
        else
            group_result(i).Waveform_peakToTrough = (waveform_t1-waveform_t_left)*interp_dt;
        end
        group_result(i).Waveform_peak = waveform_t1;
        group_result(i).Waveform_trough1 = waveform_t_left;
        group_result(i).Waveform_trough2 = waveform_t_right;
        
        % Abnormal waveform (need to be flipped)
        if group_result(i).Waveform_peakToTrough > 0.8  % 0.8 ms was set manually by plotting the distribution
            waveform = 1-waveform;
            waveform_t1 = find(waveform == max(waveform)); % Spike peak
            waveform_t_left = find(waveform == min(waveform(1:waveform_t1)),1); % Left trough
            waveform_t_right = find(waveform == min(waveform(waveform_t1:end)),1); % Right trough
            
            group_result(i).Waveform_peakToTrough = (waveform_t_right - waveform_t1)*interp_dt;
            group_result(i).Waveform_peak = waveform_t1;
            group_result(i).Waveform_trough1 = waveform_t_left;
            group_result(i).Waveform_trough2 = waveform_t_right;
            fprintf('Waveform flipped: %s\n',group_result(i).cellID{1}{1});
        end
        
        group_result(i).Waveform_broad = group_result(i).Waveform_peakToTrough > 0.35; % Set manually by plotting the distribution
        
    else
        fprintf('No waveform: %s\n', group_result(i).cellID{1}{1});
        group_result(i).Waveform_peakToTrough = nan;
        group_result(i).Waveform_broad = nan;
    end
    
    group_result(i).Waveform = waveform;
    
    % 7) Cell location.
    
    
    % 8) for dPCA
        % ZZ @ 20230704
        dPCA_data(i).PSTH = group_result(i).mat_raw_PSTH.PSTH;
        dPCA_data(i).align_markers = group_result(i).mat_raw_PSTH.align_markers;
        dPCA_data(i).align_offsets_others = group_result(i).mat_raw_PSTH.align_offsets_others;
        dPCA_data(i).trialInfo = [group_result(i).mat_raw_PSTH.stim_type_per_trial' group_result(i).mat_raw_PSTH.heading_per_trial' ...
            group_result(i).mat_raw_PSTH.choice_per_trial group_result(i).mat_raw_PSTH.outcome_per_trial];
end

% ------ Load Gaussian velocity -----
% Measured by accelerometer at 109. Should retest it on 103. HH20150422
% I've done that, they are almost the same.
temp = load('D:\Paper_rawdata\Raw_data\CN\Gaussian_vel_real_sigma045amp12.mat');
Gauss_vel = temp.Gaussian_vel_real_sigma045amp12;
Gauss_vel(:,2) = Gauss_vel(:,2)/max(Gauss_vel(:,2));


%% Reorganize data into matrices which are easier to plot and compare
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
smooth_factor_for_divergence = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(group_result);
rate_ts = group_result(representative_cell).rate_ts;
CP_ts = group_result(representative_cell).CP_ts;
memsac_ts = group_result(representative_cell).MemSac_ts;

% Initialization

% j-sensitive
% for j = 1:2
for j = 1:3
    
    % For normalized PSTH plotting (weighted by dynamic ranges for each neuron)
    PSTH_all_Norm{j} = NaN(N,length(rate_ts{j}),6);
    PSTH_correct_angles_Norm{j} = NaN(N,length(rate_ts{j}),1+length(group_result(representative_cell).unique_heading),3);  % Two zero headings
    PSTH_correctNwrong_angles_Norm{j} = NaN(N,length(rate_ts{j}),length(group_result(representative_cell).unique_heading),3);  % One zero heading
    PSTH_outcomes_Norm{j} = NaN(N,length(rate_ts{j}),4,3);
    PSTH_wrong_angles_Norm{j} = NaN(N,length(rate_ts{j}),length(group_result(representative_cell).unique_heading)-1,3); % @HH20150523
    
    % Also pack raw data for different ways of weighted sum PSTH plotting (weighted by SVM / targeted dimensionality reduction, etc.)
    PSTH_all_raw{j} = NaN(N,length(rate_ts{j}),6);
    PSTH_correct_angles_raw{j} = NaN(N,length(rate_ts{j}),1+length(group_result(representative_cell).unique_heading),3);
    PSTH_correctNwrong_angles_raw{j} = NaN(N,length(rate_ts{j}),length(group_result(representative_cell).unique_heading),3);  % One zero heading
    PSTH_outcomes_raw{j} = NaN(N,length(rate_ts{j}),4,3);
    PSTH_wrong_angles_raw{j} = NaN(N,length(rate_ts{j}),length(group_result(representative_cell).unique_heading)-1,3); % @HH20150523
    PSTH_hard_easy_raw_cellBtB4Bk{j} = NaN(N,length(rate_ts{j}),4,3); % HH20160905 For EI calculation. Cell by Time by 4 (diff, easy) by stimtype
    
    
    CP{j} = NaN(N,length(CP_ts{j}),3);
    
    ChoiceDiv_All{j} = NaN(N,length(rate_ts{j}),3);
    ChoiceDiv_All_perm{j}.std = NaN(N,length(rate_ts{j}),3);  % HH20180608
    ChoiceDiv_All_perm{j}.p = NaN(N,length(rate_ts{j}),3);
    
    ChoiceDiv_Easy_All{j} = NaN(N,length(rate_ts{j}),3);
    ChoiceDiv_Difficult_All{j} = NaN(N,length(rate_ts{j}),3);
    ChoiceDiv_EasyMinusDifficult_All{j} = NaN(N,length(rate_ts{j}),3);
    
    % Added by ZZ @20230630
    ChoiceDiv_Correct_All{j} = NaN(N, length(rate_ts{j}),3);
    ChoiceDiv_Error_All{j} = NaN(N, length(rate_ts{j}),3);
    ChoiceDiv_CorrectMinusError_All{j} = NaN(N, length(rate_ts{j}),3);
    
    ChoiceDiv_ModDiffer{j} = NaN(N,length(rate_ts{j}),3);  % 3-1, 3-2, 1-2
    ModDiv_All{j} = NaN(N,length(rate_ts{j}),3); % HH20140415   % 2-1, 3-1, 3-2
end

% j-insensitive
group_PREF_target_location = NaN(N,1);
group_PREF_target_location_notebook = NaN(N,1); % From notebook @HH20160906
group_MemSac_DDI = NaN(N,6);
group_MemSac_ps = NaN(N,6);
group_MemSac_AI = NaN(N,6);
group_TwoPtMemSac_DDI = NaN(N,6);
group_TwoPtMemSac_ps = NaN(N,6);
group_TwoPtMemSac_AI = NaN(N,6);
group_MemSac_PREF_Null_DI = NaN(N,6);
group_MemSac_actual_DI = NaN(N,1); % DI of actuall target locations (from 8pt Mem). @HH20150524
group_MemSac_actual_DI_2pt = NaN(N,1); % DI (from 2pt Mem)
group_MemSac_PREFmNULL_LeftRight_PSTH = NaN(N,length(group_result(representative_cell).MemSac_ts));
group_MemSac_PSTH_AngDiff = NaN(N,6);

group_ChoicePreference_pvalue = reshape([group_result(:).ChoicePreference_pvalue]',3,[],3);  % I updated ChoicePref time window. HH20160918
group_ChoicePreference = reshape([group_result(:).ChoicePreference]',3,[],3);  % Stim, Cell No, Pre/Post. Updated time win HH20160918
% group_ChoicePrefernce_pvalue = squeeze(group_ChoicePrefernce_pvalue(3,:,:))'; % HH20160918. 3: stim-on to stim-off

% Reorganize data
for i = 1:N
    
    %============================== NOTE ====================================
    % Changed by ZZ from normalizing each section separately to
    % normalizing the whole PSTH use the same gain !!!!
    % ========================================================
    PSTH_all_Norm_this = cellfun(@(x) x.ys, group_result(i).mat_raw_PSTH.PSTH(:,ALL_CorrectCHOICE,1), 'UniformOutput',0);
    PSTH_all_raw_this = PSTH_all_Norm_this;
    
    % FRs whthin the whole trial time share a same gian for each cell
    % FRs are normalized to [0 1]
    offset = min(cellfun(@(x) min(min(x)), PSTH_all_Norm_this));
    gain = max(cellfun(@(x) max(max(x)), PSTH_all_Norm_this)) - offset;
    
    % j-sensitive
    %     for j = 1:2
    for j = 1:3
        
        PSTH_all_Norm_this{j} = PSTH_all_Norm_this{j} - offset;
        PSTH_all_Norm_this{j} = PSTH_all_Norm_this{j} / gain;
        
        % Save weights of classical normalization method for comparison
        % with weights from other methods (targed dim reduction, etc.)
        weights_normalized_PSTH(i) = 1/gain;
        
        for stim_type = 1:3 % Stim_type check
            
            % I have done "always output three conditions" in TEMPO_GUI for all data except CP and PSTH. @HH20150419
            k = find(stim_type == group_result(i).mat_raw_PSTH.unique_stim_type);
            if ~isempty(k)   % We have this condition
                
                % ----------- Pack PSTH_all_Norm -------------
                
                PSTH_all_Norm{j}(i,:,k*2-1) = PSTH_all_Norm_this{j}(k*2-1,:);
                PSTH_all_Norm{j}(i,:,k*2) = PSTH_all_Norm_this{j}(k*2,:);
                PSTH_all_raw{j}(i,:,k*2-1) = PSTH_all_raw_this{j}(k*2-1,:);
                PSTH_all_raw{j}(i,:,k*2) = PSTH_all_raw_this{j}(k*2,:);
                
                % ---------- Normalize and pack PSTH_angles_Norm ---------
                PSTH_correct_angles_norm_this = group_result(i).mat_raw_PSTH.PSTH{j,CORRECT_ANGLE,k}.ys;
                PSTH_correct_angles_raw_this = PSTH_correct_angles_norm_this;
                PSTH_correct_angles_norm_this = PSTH_correct_angles_norm_this - offset;
                PSTH_correct_angles_norm_this = PSTH_correct_angles_norm_this / gain;
                
                PSTH_correctNwrong_angles_raw_this =  group_result(i).mat_raw_PSTH.PSTH{j,CORRECTNWRONG_ANGLE,k}.ys;
                PSTH_correctNwrong_angles_norm_this =  PSTH_correctNwrong_angles_raw_this - offset;
                PSTH_correctNwrong_angles_norm_this = PSTH_correctNwrong_angles_norm_this / gain;
                
                if size(PSTH_correct_angles_norm_this,1) == size(PSTH_correct_angles_Norm{j},3)
                    PSTH_correct_angles_Norm{j}(i,:,:,k) = PSTH_correct_angles_norm_this';
                    PSTH_correct_angles_raw{j}(i,:,:,k) = PSTH_correct_angles_raw_this';
                    
                    PSTH_correctNwrong_angles_Norm{j}(i,:,:,k) = PSTH_correctNwrong_angles_norm_this';
                    PSTH_correctNwrong_angles_raw{j}(i,:,:,k) = PSTH_correctNwrong_angles_raw_this';
                    
                elseif size(PSTH_correct_angles_norm_this,1) == size(PSTH_correct_angles_Norm{j},3) - 2 % Without zero heading
                    PSTH_correct_angles_Norm{j}(i,:,3:end,k) = PSTH_correct_angles_norm_this';
                    PSTH_correct_angles_raw{j}(i,:,3:end,k) = PSTH_correct_angles_raw_this';
                    
                    PSTH_correctNwrong_angles_Norm{j}(i,:,[1:fix(end/2) fix(end/2)+2:end],k) = PSTH_correctNwrong_angles_norm_this';
                    PSTH_correctNwrong_angles_raw{j}(i,:,[1:fix(end/2) fix(end/2)+2:end],k) = PSTH_correctNwrong_angles_raw_this';
                else
                    disp('No match PSTH_angles_norm...');
                end
                
                if group_result(i).mat_raw_PSTH.PREF == LEFT  % Note that PSTH{6} is grouped by raw heading, not flip to "RIGHT = PREF"
                    PSTH_correctNwrong_angles_Norm{j}(i,:,:,k) = fliplr(squeeze(PSTH_correctNwrong_angles_Norm{j}(i,:,:,k)));
                    PSTH_correctNwrong_angles_raw{j}(i,:,:,k) = fliplr(squeeze(PSTH_correctNwrong_angles_raw{j}(i,:,:,k)));
                end
                
                PSTH_outcome_norm_this = group_result(i).mat_raw_PSTH.PSTH{j,OUTCOME,k}.ys;
                PSTH_outcome_raw_this = PSTH_outcome_norm_this;
                
                PSTH_outcome_norm_this = PSTH_outcome_norm_this - offset;
                PSTH_outcome_norm_this = PSTH_outcome_norm_this / gain;
                
                PSTH_outcomes_Norm{j}(i,:,:,k) = PSTH_outcome_norm_this';
                PSTH_outcomes_raw{j}(i,:,:,k) = PSTH_outcome_raw_this';
                
                PSTH_wrong_angles_norm_this = group_result(i).mat_raw_PSTH.PSTH{j,WRONG_ANGLE,k}.ys;
                order = nan(size(group_result(i).mat_raw_PSTH.sort_info{WRONG_ANGLE}{2}{2,3}));
                
                if group_result(i).mat_raw_PSTH.PREF == 1 % If PREF is leftward, then the positive headings go first (because they are WRONG trials)
                    order(1:2:end) = length(order)/2+1:1:length(order);
                    order(2:2:end) = length(order)/2:-1:1;
                else % vice versa
                    order(1:2:end) = length(order)/2:-1:1;
                    order(2:2:end) = length(order)/2+1:1:length(order);
                end
                
                PSTH_wrong_angles_norm_this = PSTH_wrong_angles_norm_this(order,:);
                PSTH_wrong_angles_raw_this = PSTH_wrong_angles_norm_this;
                
                PSTH_wrong_angles_norm_this = PSTH_wrong_angles_norm_this - offset;
                PSTH_wrong_angles_norm_this = PSTH_wrong_angles_norm_this / gain;
                
                PSTH_wrong_angles_Norm{j}(i,:,:,k) = PSTH_wrong_angles_norm_this';
                PSTH_wrong_angles_raw{j}(i,:,:,k) = PSTH_wrong_angles_raw_this';
                
                % HH20160905 EI
                PSTH_hard_easy_raw_cellBtB4Bk{j}(i,:,:,k) = group_result(i).mat_raw_PSTH.PSTH{j,CHOICE_DIFFICULT,k}.ys';
            end
            
            % Stim type check already done
            ChoiceDiv_All{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_ALL{j}(stim_type,:),smooth_factor_for_divergence);
            ChoiceDiv_All_perm{j}.std(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_ALL_perm{j}.std(stim_type,:),smooth_factor_for_divergence);
            ChoiceDiv_All_perm{j}.p(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_ALL_perm{j}.p(stim_type,:),smooth_factor_for_divergence);
            
            % Difficult and easy trials
            ChoiceDiv_Easy_All{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_Easy{j}(stim_type,:),smooth_factor_for_divergence);
            ChoiceDiv_Difficult_All{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_Difficult{j}(stim_type,:),smooth_factor_for_divergence);
            ChoiceDiv_EasyMinusDifficult_All{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_Easy{j}(stim_type,:) ...
                - group_result(i).mat_raw_PSTH.ChoiceDivergence_Difficult{j}(stim_type,:),smooth_factor_for_divergence);
            
            % Correct and error trials
            ChoiceDiv_Correct_All{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_Correct{j}(stim_type,:), smooth_factor_for_divergence);
            ChoiceDiv_Error_All{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_Error{j}(stim_type,:), smooth_factor_for_divergence);
            ChoiceDiv_CorrectMinusError_All{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ChoiceDivergence_Correct{j}(stim_type,:) ...
                - group_result(i).mat_raw_PSTH.ChoiceDivergence_Error{j}(stim_type,:), smooth_factor_for_divergence);
            
            CP{j}(i,:,stim_type) = group_result(i).CP{j}(stim_type,:);
            CP_p{j}(i,:,stim_type) = group_result(i).CP_p_value{j}(stim_type,:);
            
            ModDiv_All{j}(i,:,stim_type) = smooth(group_result(i).mat_raw_PSTH.ModalityDivergence{j}(stim_type,:),smooth_factor_for_divergence);
            
        end
        
        ChoiceDiv_ModDiffer{j}(i,:,1) = smooth(ChoiceDiv_All{j}(i,:,3) - ChoiceDiv_All{j}(i,:,1),smooth_factor_for_divergence);
        ChoiceDiv_ModDiffer{j}(i,:,2) = smooth(ChoiceDiv_All{j}(i,:,3) - ChoiceDiv_All{j}(i,:,2),smooth_factor_for_divergence);
        ChoiceDiv_ModDiffer{j}(i,:,3) = smooth(ChoiceDiv_All{j}(i,:,1) - ChoiceDiv_All{j}(i,:,2),smooth_factor_for_divergence);
        % Different with HH, Modality preference I calculate using all correct and wrong
        % trials. And pairs are 12, 13, 23
        % Commented by ZZ 20210506
    end
    % Aligned Choice Preference 
    ChoicePref_Difficult_All(i,:,:) = squeeze(group_result(i).mat_raw_PSTH.ChoicePref_difficeasy(:,:,1)) * sign(group_result(i).if_contralateral - 0.5);
    ChoicePref_Easy_All(i,:,:) = squeeze(group_result(i).mat_raw_PSTH.ChoicePref_difficeasy(:,:,2))* sign(group_result(i).if_contralateral - 0.5);
    ChoicePref_Correct_All(i,:,:) = squeeze(group_result(i).mat_raw_PSTH.ChoicePref_correcterror(:,:,1))* sign(group_result(i).if_contralateral - 0.5);
    ChoicePref_Error_All(i,:,:) = squeeze(group_result(i).mat_raw_PSTH.ChoicePref_correcterror(:,:,2))* sign(group_result(i).if_contralateral - 0.5);

    % j-insensitive
        
    % --- Find the PREF target location ---
    % (1) From the eye trace. See HeadingDis_cum_PSTH_HH.m
    tmp = group_result(i).mat_raw_PSTH.PREF_target_location; % [-90,270]. There was a bug:)
    group_PREF_target_location(i) = mod(tmp,360);  % Changed to [0,360] @HH20160906
    
    % (2) Use the notebook. Note I only record the right target (not necessarily the PREF target). HH20160906
    tmp = xls_num{1}(i,header.HD_TargAng); % [-90,90]
    if group_result(i).PREF_PSTH == LEFT  % Flip 180 degree
        tmp = tmp + 180;   % [-90,270]
    end
    group_PREF_target_location_notebook(i) = mod(tmp,360); % [0,360]
    
    
    % Mem-sac stuff (Note some are NaNs)
    if ~isempty(group_result(i).MemSac_DDI)
        % -------- Global memsac indicator -------
        group_MemSac_DDI(i,:) = group_result(i).MemSac_DDI;
        group_MemSac_ps(i,:) = group_result(i).MemSac_p;
        group_MemSac_AI(i,:) = group_result(i).MemSac_AI;
        
        % -------- Local memsac indicator --------
        % 1. Left and Right (Pref/null of HD task)
        group_MemSac_PREFmNULL_LeftRight_PSTH(i,:) = (group_result(i).MemSac_PSTH(:,1) - group_result(i).MemSac_PSTH(:,5)) *  sign((group_result(i).PREF_PSTH == 2)-0.5);  % (Right - Left)* Right is pref
        group_MemSac_PSTH_AngDiff(i,:) = mod(group_result(i).MemSac_vectSum - ((group_result(i).PREF_PSTH == 2)*0 + (group_result(i).PREF_PSTH == 1)*180),360);
        
        % 2. Mem's Pref_null DI. @HH20150524
        group_MemSac_PREF_Null_DI(i,:) = group_result(i).mat_raw_MemSac.PREF_NULL_DI;
        
        % 3. Actual DI of PREf/NULL of HD task.  @HH20150524
        
        % Interpolate original Memsac traces into higher spatial resolution
        
        % Save data for individual cell plotting
        group_result(i).MemSac_interp_locations = group_result(i).mat_raw_MemSac.MemSac_interp_locations;
        group_result(i).MemSac_interp_PSTH = group_result(i).mat_raw_MemSac.MemSac_interp_PSTH{3};
        group_result(i).MemSac_interp_PSTH_alignVisON = group_result(i).mat_raw_MemSac.MemSac_interp_PSTH{1}; % Added align to Vis Off
        
        
        % --- Calculate Memsac DI using user defined period ----
        MemSac_actual_DI_period = [3];  % [3,4] % I use memory and pre period to calculate the actual DI
        
        MemSac_temporal_Slice = group_result(i).mat_raw_MemSac.temporal_Slice;
        MemSac_align_offsets = group_result(i).mat_raw_MemSac.align_offsets;
        MemSac_align_markers = group_result(i).mat_raw_MemSac.align_markers;
        MemSac_ts = group_result(i).MemSac_ts;
        MemSac_actual_DI_time_ind = false(size(MemSac_ts));
        
        for p_ind = 1:length(MemSac_actual_DI_period)
            
            ppp = MemSac_actual_DI_period(p_ind);
            
            if MemSac_temporal_Slice{ppp,3} == 7  % Precise windows (because MemSac_ts has been aligned to 7 (saccade onset))
                add_ind = MemSac_temporal_Slice{ppp,1} <= MemSac_ts & MemSac_ts <= MemSac_temporal_Slice{ppp,2};
            else  % Windows that is not so precise
                % Mean shift
                meanShift = mean(MemSac_align_offsets(:,MemSac_align_markers == MemSac_temporal_Slice{ppp,3})-MemSac_align_offsets(:,MemSac_align_markers==7),1);
                add_ind = meanShift + MemSac_temporal_Slice{ppp,1} <= MemSac_ts & MemSac_ts <= meanShift + MemSac_temporal_Slice{ppp,2};
            end
            
            MemSac_actual_DI_time_ind = MemSac_actual_DI_time_ind | add_ind ;
        end
        
        % Merge the gap
        MemSac_actual_DI_time_ind(find(MemSac_actual_DI_time_ind,1):find(MemSac_actual_DI_time_ind,1,'last')) = true;
        
        % Use the interpolated data from 8p Memsac
        % Find actual DI for this cell
        [~,pref] = min(abs(group_PREF_target_location(i) - group_result(i).MemSac_interp_locations));
        pref = mod(pref-1,length(group_result(i).MemSac_interp_locations)-1)+1;
        null = mod(pref + (length(group_result(i).MemSac_interp_locations)-1)/2 -1, length(group_result(i).MemSac_interp_locations)-1)+1; % The opposite position
        
        MemSac_PREF_mean = mean(group_result(i).MemSac_interp_PSTH(MemSac_actual_DI_time_ind,pref),1);
        MemSac_NULL_mean = mean(group_result(i).MemSac_interp_PSTH(MemSac_actual_DI_time_ind,null),1);
        
        group_MemSac_actual_DI(i) = (MemSac_PREF_mean - MemSac_NULL_mean) / (MemSac_PREF_mean + MemSac_NULL_mean);
        
    end
    
    group_MemSac_PSTH_AngDiff(group_MemSac_PSTH_AngDiff > 180) = 360 - group_MemSac_PSTH_AngDiff(group_MemSac_PSTH_AngDiff > 180); % Ang diffs are less than 180
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------- Choose Memsac indicator (affect all Memsac measures) -------------
% 1 = Background, 2 = LS, 3 = Mem, 4 = Pre, 5 = Co, 6 = Post

MemSac_indicator = mean(group_MemSac_DDI(:,[3]),2); % Global DDI
MemSac_indicator_p = group_MemSac_ps(:,3);
MemSac_indicator_txt = 'MemSac\_DDI ([3])';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Psychometric parameters
Psy_all = [group_result.Psy_para];
Psy_thres = reshape(Psy_all(1,:),3,[])';

bayes_pred = @(x,y) sqrt(1./(1./x.^2 + 1./y.^2));
Psy_thres_pred = bayes_pred(Psy_thres(:,1),Psy_thres(:,2));
Psy_pred_ratio = Psy_thres(:,3)./Psy_thres_pred;
Psy_pred_ratio_vestibular_visual = Psy_thres(:,1:2)./repmat(Psy_thres_pred,1,2);


select_all = []; select_sus = []; select_bottom_line = []; select_tcells = []; select_no_tcells = []; select_bottom_line_all_monkey = [];
select_psy_good = []; select_psy_bad = []; find_bottom_line = [];

Choice_pref_all = []; Choice_pref_p_value_all = []; Modality_pref_all = []; Modality_pref_p_value_all = [];
select_cpref_mpref = [];

t_criterion_txt = [];

group_position = [];  % All cells
Position_all = []; % Could be one monkey
cell_position();
cell_selection();
% Commented by ZZ 20201207

%% Cell Selection and Cell Counter
    function cell_selection(t_cell_selection_num)  % Cell Selection and Cell Counter
        
        if nargin < 1
            t_cell_selection_num = 6; % Default
        end
        
        % ---------- ZZ @20231012 ------------
        % Excluding cells with SNR < 6
        SNR_ITI = xls_txt{1}(:, header.SNR_ITI);    % Two values, the first one is SNR and the latter is median(ITI)/mean(ITI)
        space_ind = cellfun(@(x) strfind(x, ' '), SNR_ITI, 'uniformoutput', 0);
        % SNR
        SNR = cell2mat(cellfun(@(x, y) str2num(x(y(1)+1:y(2)-1)), SNR_ITI, space_ind, 'uniformoutput', 0));
        snr_large_ind = SNR > 6; 
        
        % Try to sort cells into different subgroups (MSN, TAN, FSN)
        % median(ISI) / mean(ISI)
        ISI = cell2mat(cellfun(@(x, y) str2num(x(y(2)+1:end)), SNR_ITI, space_ind, 'uniformoutput', 0));
        % waveform related
        Peak2Trough = ([group_result.Waveform_peakToTrough])';
        % Spontaneous FR
        % Defined as the mean FR during 600~800ms after Feedback (10ms window, 5ms step)
        time_start = 600:5:790;
        for i = 1:length(group_result)
            for n = 1:length(time_start)
                Spon_FR_all(i,n) = sum(mean(group_result(i).mat_raw_PSTH.spike_aligned{1,3}(:,time_start(n):time_start(n)+9))) / 10*1000;
            end
        end
        Spon_FR = mean(Spon_FR_all,2);
        
        % 3 variables for classification 
        var_class = [log(ISI) Peak2Trough Spon_FR];   
        
        % I use AO for recording latter, this may affect Peak2Trough
        ced_data = [group_result.cell_num]' < 200;
        ao_data = [group_result.cell_num]' >200;
        
        [~, ~, ~, ced_outliner] = robustcov(var_class(ced_data,:), 'outlierfraction', 0.05); 
        [~, ~, ~, ao_outliner] = robustcov(var_class(ao_data,:), 'outlierfraction', 0.2); 
        outlier = [ced_outliner; ao_outliner];
        
        % Plotting
        temp_ind = logical(zeros(sum(ced_data),1)); 
%         set(figure(1012), 'name', 'Classification of Cells');clf; hold on; 
%         plot3(var_class(ced_data,1), var_class(ced_data,2), var_class(ced_data,3), 'ko', 'markersize', 10);
%         plot3(var_class(ao_data,1), var_class(ao_data,2), var_class(ao_data,3), 'k^', 'markersize', 10);
%         
%         plot3(var_class(ced_outliner,1), var_class(ced_outliner,2), var_class(ced_outliner,3), 'ko', 'markerfacecolor','k', 'markersize', 10);
%         plot3(var_class([temp_ind;ao_outliner],1), var_class([temp_ind;ao_outliner],2), var_class([temp_ind;ao_outliner],3),...
%             'k^',  'markerfacecolor','k','markersize', 10);
%         
%         xlabel('ITI'); ylabel('Waveform'); zlabel('Firing rates'); 
%         title([num2str(length(var_class)) ' cells with ' num2str(sum(outlier)) ' outliers']); 
%         SetFigure();
        
        
        % --------  @ HH20150413 --------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Limitations on petition number & Target first
        select_all_all_monkey = ([group_result.repN]' >= 8) & ([group_result.length_unique_stim_type]' == 3 ) & (xls_num{1}(:,header.HD_TargFirst)~=0)& snr_large_ind & outlier~=1;  % Added by ZZ @ 20231012
        % + SUs
        select_sus_all_monkey = select_all_all_monkey & (xls_num{1}(:,header.Units_RealSU) > 0) & (xls_num{1}(:,header.Chan1) < 20);
        
        % Bottom line for most figures
        select_bottom_line_all_monkey = select_sus_all_monkey;
        
        % + T(ypical) Cells
        t_cell_selection_criteria = ...
            {  %                   Logic                                   Notes
            % Bottom-line
            'Bottom-line (all)', ones(size(select_sus_all_monkey)); % Just all bottom-line cells
            % Mem-sac based
            'Manually assigned mem-sac (original)',    xls_num{1}(:,header.HD_MemSac) >= 0.8;
            'Memory p value',   group_MemSac_ps(:,3) < 0.05 | group_TwoPtMemSac_ps(:,3) < 0.05;
            'Memory p value & actual DDI', (group_MemSac_ps(:,3)<0.05 & abs(group_MemSac_actual_DI) >= 0.3) | (group_TwoPtMemSac_ps(:,3)<0.05 & abs(group_TwoPtMemSac_DDI(:,3)) >=0.3);
            'Global DDI of mem-sac',   mean(group_MemSac_DDI(:,[3 4]),2) >= 0.5;
            % HD based
            'ChoicePreference p value (any)', any(group_ChoicePreference_pvalue(:,:,1) < 0.01,1)';  % 3, stim-on to stim-off
            'ChoicePreference p value (all)', all(group_ChoicePreference_pvalue(:,:,1) < 0.01,1)';
            'ChoicePreference p value (vest)', (group_ChoicePreference_pvalue(1,:,1) < 0.01)';  % 3, stim-on to stim-off
            'ChoicePreference p value (vis)', (group_ChoicePreference_pvalue(2,:,1) < 0.01)';  % 3, stim-on to stim-off
            'ChoicePreference p value (comb)', (group_ChoicePreference_pvalue(3,:,1) < 0.01)';  % 3, stim-on to stim-off
            %             'ModalityPreference p value (any)', any(group_ModalityPreference_pvalue(:,:,1) < 0.01,1)
            % For debugging
            'HD_dt patched', [group_result(:).HD_dt_patched]';
            'non HD_dt patched', ~[group_result(:).HD_dt_patched]';
            'Combined > Max (Vis, Vest) and Comb p < 0.01',(group_ChoicePreference_pvalue(3,:,1) < 0.01)' & (abs(group_ChoicePreference(3,:,1)) > abs(group_ChoicePreference(1,:,1)))' & (abs(group_ChoicePreference(3,:,1)) > abs(group_ChoicePreference(2,:,1)))';
            };
        
        select_tcells_all_monkey = select_sus_all_monkey & t_cell_selection_criteria{t_cell_selection_num,2};
        select_no_tcells_all_monkey = select_sus_all_monkey & ~ t_cell_selection_criteria{t_cell_selection_num,2};
        t_criterion_txt = t_cell_selection_criteria{t_cell_selection_num,1};
        
        % ---------
        selected_t_vest = select_sus_all_monkey & (group_ChoicePreference_pvalue(1,:,1) < 0.01)';
        selected_t_vis = select_sus_all_monkey & (group_ChoicePreference_pvalue(2,:,1) < 0.01)';
        selected_t_comb = select_sus_all_monkey & (group_ChoicePreference_pvalue(3,:,1) < 0.01)';
        
        % Psychometric good (prediction ratio <= 1.3)
        select_psy_good = Psy_pred_ratio <= 1.3;
        select_psy_bad = Psy_pred_ratio > 1.3;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % -------- Count cell numbers for each monkey. HH20150723 --------
        n_monkey = length(monkey_included_for_loading);
        cell_nums = zeros(n_monkey + 1,3); % All units / SUs / T Cells
        for mm = 1:n_monkey
            select_monkey{mm} = xls_num{1}(:,header.Monkey) == monkey_included_for_loading(mm);
            cell_nums(mm,1) = sum(select_all_all_monkey & select_monkey{mm});
            cell_nums(mm,2) = sum(select_sus_all_monkey & select_monkey{mm});
            cell_nums(mm,3) = sum(select_tcells_all_monkey & select_monkey{mm});
        end
        
        % -------- Update actual dataset for analysis. HH20150723 --------
        monkey_included_for_analysis = monkey_included_for_loading(logical([get(findall(gcbf,'tag','Duncan_data'),'value') get(findall(gcbf,'tag','Fara_data'),'value')]));
        %         monkey_included_for_analysis = monkey_included_for_loading(logical(get(findall(gcbf,'tag','Duncan_data'),'value')));
        monkey_mask_for_analysis = false(length(group_result),1);
        for mm = 1:length(monkey_included_for_analysis)
            monkey_mask_for_analysis = monkey_mask_for_analysis | (xls_num{1}(:,header.Monkey) == monkey_included_for_analysis(mm));
        end
        
        % -------- Affect all analysis below --------
        select_all = select_all_all_monkey & monkey_mask_for_analysis;
        select_sus = select_sus_all_monkey & monkey_mask_for_analysis;
        select_bottom_line = select_bottom_line_all_monkey & monkey_mask_for_analysis;   find_bottom_line = find(select_bottom_line);
        select_tcells = select_tcells_all_monkey & monkey_mask_for_analysis;
        select_no_tcells = select_no_tcells_all_monkey & monkey_mask_for_analysis;
        
        cell_nums(end,:) = [sum(select_all) sum(select_sus) sum(select_tcells)];
        
        % -------- Update cell counter ---------
        h_all = findall(gcbf,'tag','num_all_units');
        set(h_all,'string',sprintf('%7d%7d%7d\n',cell_nums'),'fontsize',13);
        h_t_criterion = findall(gcbf,'tag','t_criterion');
        set(h_t_criterion,'string',{t_cell_selection_criteria{:,1}});
        set(h_t_criterion,'value',t_cell_selection_num);
        
        % -------- Update/refresh some related datasets that are influenced by cell_selection ---------
        % For Choice and modality preference
        select_cpref_mpref = select_bottom_line;
        Choice_pref_all = reshape([group_result(select_cpref_mpref).ChoicePreference]',3,[],3);  % Stim, Cell No, Pre/Post. Updated time win HH20160918
        Choice_pref_p_value_all = reshape([group_result(select_cpref_mpref).ChoicePreference_pvalue]',3,[],3); % Updated time win. HH20160918
        Modality_pref_all = reshape([group_result(select_cpref_mpref).ModalityPreference]',3,[],3);
        Modality_pref_p_value_all = reshape([group_result(select_cpref_mpref).ModalityPreference_pvalue]',3,[],3);
        
        ChoicePref_Difficult = ChoicePref_Difficult_All(select_cpref_mpref,:,:) ; 
        ChoicePref_Easy = ChoicePref_Easy_All(select_cpref_mpref,:,:); 
        ChoicePref_Correct = ChoicePref_Correct_All(select_cpref_mpref,:,:);
        ChoicePref_Error = ChoicePref_Error_All(select_cpref_mpref,:,:);
        
        select_group_position = group_position(select_cpref_mpref,:); 
        
        % Added by ZZ @ 20230703
        for j = 1:3 
            ChoiceDiv_All_perm_select{j}.p = ChoiceDiv_All_perm{j}.p(select_cpref_mpref,:,:);
            
            ChoiceDiv_Error{j} = ChoiceDiv_Error_All{j}(select_cpref_mpref,:,:);
            ChoiceDiv_Correct{j} = ChoiceDiv_Correct_All{j}(select_cpref_mpref,:,:);
            ChoiceDiv_CorrectMinusError{j} = ChoiceDiv_CorrectMinusError_All{j}(select_cpref_mpref,:,:);
            ChoiceDiv_Difficult{j} = ChoiceDiv_Difficult_All{j}(select_cpref_mpref,:,:);
            ChoiceDiv_Easy{j} = ChoiceDiv_Easy_All{j}(select_cpref_mpref,:,:);
            ChoiceDiv_EasyMinusDifficult{j} = ChoiceDiv_EasyMinusDifficult_All{j}(select_cpref_mpref,:,:);
        end
        
        select_for_SVM = select_bottom_line;
        select_for_PCA_B = select_bottom_line;
        
        PCA_A = [];  % Reset PCA_A
        PCA_B_projPC = []; % Reset PCA_B
        thres_choice = []; % Reset SVM training
        weights_PCA_B_PC = []; weights_svm_choice_mean = []; weights_TDR_PCA_SVM_mean = []; % Reset TDR
        
    end
% Commented by ZZ 20201207


%% Load DrawMapping data to get the cell location. HH20161024
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
            
        [15,1],[31,0];       % Duncan_L    % Changed from AP0 to AC0
        [15,2],[31,0];       % Duncan_R
        [13,1],[19,0];       % Fara_L
        [13,2],[19,0];       % Fara_R
        %             [15,1],[55,0]
        %             [15,2],[55,0]
        
        };
    
    drawmapping_data{1,3} = {
        %{[Session(s)], [LocX(Posterior) LoxY(Lateral)], [GuideTube(cm) Offset(cm)], [AreaType, Begin(100um), End(100um); ...] , electrode retrieval}
        {10,[30,9],[1.4,0.0],[FC 20 48; CD 148 195]}
        {11,[30,7], [1.4,-0.2], [FC 14 36; CC 76 101; CD 144 184; GP 210 230]}    % 从guidetube内距离尖端2mm处出发
        {12,[30,11],[1.8,0.0],[FC 6 27; CD 160 191; GP 200 240; GM 260 270]}
        {13, [30,5], [1.75,-0.1], [FC 23 36; CC 58 80.5; CC 86 104; CD 138 195; GP 210 230]}  %约210进入GP 但是不知道多少出230，这里为了不报错随便写了一个值
        %             {14, [30,7], [1.5,0.0], [FC 14 36; CC 76 101; CD 144 184; GP 210 230]}  % Right hemisphere
        {15, [30,7], [1.75,-0.1], [FC 23 36; CC 58 80.5; CC 86 104; CD 138 195; GP 210 230]}
        {16, [30,11],[1.75 -0.1], [FC 21 48; CD 150 163; Pu 195 207; GP 219 239]}
        {17, [30,6], [1.5 0.0], [FC 15 40; CC 73 97; CC 92 108; CD 149 197; GP 230 250]}
        {18, [30,5], [1.75,-0.1], [FC 23 36; CC 58 80.5; CC 86 104; CD 138 195; GP 210 230]}
        {19, [30,8], [1.6,0.0], [FC 9 34; CD 136 167]} % 只走到了167
        {20, [29,7], [1.6,0.0], [CC 76 105; CD 144 172]} % 只到了172
        {21, [32,7], [1.75,0.0],[FC 9 23; CC 61 91; CD 127 149]} % only to 149
        {22, [35,8], [1.75,0.0],[FC 4 30; CC 66 83; CD 135 158]} % 158 out of CD
        {23, [38,9], [1.75,0.0],[FC 7 30; CD 128 141]}
        {24, [40,9], [1.75,0.0],[FC 7 29; CD 136 148]}
        {25, [30,9], [1.75,0.0],[FC 10 28; CD 119 150]}  % only to 144
        {26, [28,6], [1.75,-0.1],[FC 13 31; CC 57 90; CD 128 170; GP 209 229]}
        {27, [26,6], [1.6,0.0], [FC 10 31; CC 63.5 89; CD 142 168]}    % only to 168
        {28, [22,6], [1.6,0.0], [FC 24 48; CC 73.7 111; CD 156 220]}   % CD= CD+NAc only to 220
        {29, [20,6], [1.75,0.0],[FC 3 25; CC 52 84; CD 132 193]} % only to 193
        {30, [31,6], [2.0,0.0], [CC 10 28; CC 33.5 48; CD 93 109]} % only to 109
        {31, [31,7], [1.75,0.0],[FC 3 31; CC 47 64.5; CC 70 86; CD 131 168]} % only to 168
        {32, [31,8], [1.75,0.0],[CC 78.5 88; CD 136 178]}
        {33, [31,9], [1.6, 0.0],[FC 16 38; CD 139 155]} % only to 155
        {34, [30,10],[1.6, 0.0],[FC 21 51; CD 140 166]}
        {35, [20,7], [1.6, 0.0],[FC 15 43; CD 151 164]} % only to 153
        {36, [20,8], [1.6, 0.0],[FC 21 57; CD 150 160]} % only to 160
        {37, [20,9], [1.6, 0.0],[FC 15 44; CD 141 176]} % only to 176
        {38, [20,10],[1.6, 0.0],[FC 13 49; FC 78 93; CD 149 167]} % only to 167
        {39, [19,7], [1.6, 0.0],[FC 14 50; CC 78 101;LV 149 154; CD 155 181]} % only to 181
        {40, [19,8], [1.6, 0.0],[FC 11 45; CD 147 167]} % only to 167
        {41, [19,9], [1.6, 0.0],[FC 18 45; FC 96 97; CD 145 179]} % only to 179
        {42, [19,10],[1.6, 0.0],[FC 21 48; CD 146 185; GM 208 228]} % only to 208
        {43, [18,7], [1.6, 0.0],[FC 11 42; CC 79 95; CD 157 205; GM 217 237]} % only to 225
        {44, [18,8], [1.6, 0.0],[FC 10 34; GM 63 95; LV 142 155; CD 159 200]} % only to 200
        {45, [18,9], [1.6, 0.0],[FC 20 63; CD 167 196]} % only to 196
        {46, [18,10],[1.6, 0.0],[FC 9 40; GM 47 71; CD 145 169]} % only to 169
        {47, [17,7], [1.75,0.0],[FC 3 25; CC 58 82; LV 135 140; CD 140 152]} % only to 152
        {48, [17,8], [1.75,0.0],[FC 7 29; CC 50 93; LV 143 160; CD 162 171; GM 211 223]}
        {49, [17,9], [1.75,0.0],[FC 3 26; CD 133 154]} % only to 154
        {50, [17,10],[1.75,0.0],[FC 7 33; CD 141 153]} % only to 153
        {51, [17,11],[1.8, 0.0],[FC 1 25; CD 129 135]} % only to 135
        {52, [17,6], [1.75,0.0],[FC 10 33; CC 53 73; CC 80 93; LV 142 158; CD 158 174; GM 197 199]} % only to 199
        {53, [18,6], [1.75,0.0],[FC 4 28; CC 58 86; LV 129 134; CD 135 150]} % only to 150
        {54, [17,12],[1.8, 0.0],[FC 5 77]}
        {54, [17,11],[1.7,0.0], [FC 22 53; GM 206 220]}
        {55, [15,6], [1.8,0.0], [FC 11 35; CC 55 72; CC 77 93; LV 141 150; CD 155 221]}
        {56, [15,7], [1.75,0.0], [FC 21 45; CC 66 99; CD 160 178; GM 192 212]}
        {57, [20,5], [1.75,0.0], [FC 14 35; CC 57 76; CC 83 94; CD 137 211]} % VM in 23 but not in near holes; CD = CD + NAc; only to 211
        {58, [20,4], [1.75,0.0], [FC 15 42; CC 60 81; CC 85 100; LV 137 148; CD 148 180; LV 183 190; NAc 190 231]} % no VM; only to 231; NAc ?
        {58, [20,11],[1.75,0.0], [FC 22 98; CD 145 160; Pu 176 196]} % only to 196
        {59, [20,3], [1.75,0.0], [FC 23 50; CC 63 87; CC 96 110; LV 149 155; CD 155 177; LV 177 194; NAc 194 212]} % only to 212
        {60, [19,5], [1.8, 0.0], [FC 23 43; CC 60 76; CC 84 98; LV 133 145; CD 146 187]} % CD = CD + NAc; only to 187
        {61, [19,11],[1.8, 0.0], [GM 42 65; GM 83 98; CD 147 174; Pu 187 206]} % only to 206
        {62, [19,6], [1.6, 0.0], [CC 68 83; CC 107 123; CD 181 201]} % only to 181
        {62, [20,7], [1.6, 0.0], [CC 56 75; CD 155 199]} % only to 199
        {63, [21,7], [1.6, 0.0], [CC 59 79; CC 101 117; CD 160 187]} % only to 187
        {64, [21,9], [1.6, 0.0], [CC 71 87; CD 160 176]} % only to 176
        {65, [21,8], [1.6, 0.0], [CC 75 90; CD 167 193]} % only to 193
        {66, [21,6], [1.6, 0.0], [CC 72 99; CC 108 129; CD 164 210]}
        {67, [21,10],[1.6, 0.0], [GM 81 101;CD 167 195]} % only to 195
        {68, [25,7], [1.6, 0.0], [CC 64 101; CD 163 191]} % only to 191
        {69, [25,9], [1.6, 0.0], [GM 74 98; CD 161 183]} % only to 183
        {70, [27,7], [1.6, 0.0], [CC 76 108; CD 186 211]} % only to 211   I'm not sure whether the mapping is wrong
        {71, [28,7], [1.6, 0.0], [CC 84 112; CD 188 204]} % only to 204
        {72, [25,8], [1.6, 0.0], [CC 84 110; CD 171 197]} % only to 197
        {73, [24,8], [1.6, 0.0], [CC 88 108; CD 167 185]} % only to 185
        {74, [26,7], [1.6, 0.0], [CC 83 110; CD 174 200]} % only to 200
        {75, [24,7], [1.6, -0.1],[CC 101 121; CD 180 195]} % only to 195
        {76, [27,8], [1.6, -0.1],[CC 82 108; CD 185 207]} % only to 207
        {77, [35,7], [1.7, 0.0], [CD 140 168]} % only to 168
        {78, [37,8], [1.7, -0.1],[CC 83 108; CD 150 180]} % only to 180
        {79, [23,7], [1.6, 0.0], [CC 56 78; CD 156 186]} % only to 186
        {80, [33,7], [1.8, 0.0], [FC 13 39; CC 70 104; CD 138 162]}
        {81, [36,8], [1.8, 0.0], [FC 18 48; CC 82 106; CD 147 177; GP 183 203]} % only to 183
        {82, [36,7], [1.8, 0.0], [FC 17 46; CC 77 108; CD 147 174; GP 184 204]} % only to 187
        {83, [36,6], [1.8, 0.0], [FC 13 40; CC 69 87; CC 94 111; CD 148 173; GP 179 228]} % till 228
        {84, [36,9], [1.9, 0.0], [FC 9 36; CD 139 151]} %only to 151
        {85, [35,7], [1.8, 0.0], [FC 16 43; CC 75 90; CC 96 110; CD 148 180; GP 186 206]}
        {86, [32,8], [2.05,0.0], [FC 8 26; CD 105 117]} % only to 117
        {87, [32,9], [2.15,0.0], [CD 103 115]}  %only to 115
        {88, [32,10],[2.05,0.0], [CD 109 133]} % only to 133
        {89, [33,8], [2.15,0.0], [CD 95 127]} % only to 127
        {90, [33,10],[2.05,0.0], [CD 109 124]} % only to 124
        {91, [35,10],[2.05,0.0], [CD 113 127]} % only to 127
        {92, [36,10],[2.05,0.0], [CD 115 140; Pu 149 169]} % 149 is around the next GM
        {93, [34,8], [2.05,0.0], [CD 113 135; Pu 143 163]} % only to 148
        {94, [34,10],[2.05,0.0], [CD 110 128]} % only to 128
        {95, [22,9], [2.15,0.0], [CD 103 132]} % only to 132
        {96, [22,10],[2.05,0.0], [CD 119 145]} % only to 145
        {97, [22, 8], [2.15, 0.0], [CD 102 142]} % only to 142
        {98, [23,11],[2.05, 0.0], [CD 117 142]} % only to 142
        {99, [24,10],[2.05, 0.0], [CD 104 145; GM 171 180]} % only to 177
        {100, [24,11],[2.05,0.0],[GM 59 79; CD 116 121]} % only to 121
        {101, [25,11],[2.05,0.0],[CD 112 122]} % only to 112
        {102, [26,11],[2.05,0.0],[CD 112 142; GM 159 169]} % only to 166
        {103,[26, 10],[2.05,0.0],[GM 12 33; GM 65 88; CD 121 135]} % only to 135
        {104,[27, 11],[2.05,0.0],[GM 5 33; CD 113 132]} % only to 132
        {105,[28, 11],[2.05,0.0],[GM 11 40; CD 114 138; GM 151 161]} % only to 152
        {117,[55, 16],[2.05,0.0],[GM 16 40; GM 70 90]}  % only to 90, modest VM in 7988
        {118,[53, 18],[2.05,0.0],[GM 23 49; GM 57 87]} % 10382 colse to MST
        {120,[55, 18],[2.15,0.0],[GM 10 42]}
        {123, [53, 18],[2.25,0.0],[GM 38 56]}  % only to 5080
        {124, [53, 19],[2.25,0.0], [GM 17 43]}  % only to 4227
        {125, [53, 19],[2.25,0.0], [GM 17 57; GM 91 95]}  %only to 9514
        {126, [53, 16],[2.15,0.0], [GM 14 52; GM 71 86]}
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
        {277, [15, 8], [2.4, 0.0], [GM 8 22; CD 124  154]}  % only to 154
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
        {323, [10, 5], [2.4, 0.0], [GM 20 41; CC 75 102; CD 151 180]} % only to 175
        {324, [9, 8], [2.4, 0.0], [GM 22 39; GM 89 107; CD 155 180]}  % only to 178
        };
    
    %             }';
    
    
    % === Get cell positions ===
    for i = 1:length(group_result)
        % Decode cell position
        id = sscanf(group_result(i).cellID{1}{1},'%gm%gs%gh%gx%gy%gd%g');
        this_monkey_hemi = id([2 4])';
        this_session = id(3);
        this_pos_raw = id(5:7)';
        
        % Find entry in drawmapping_data
        found = 0;
        for mmhh = 1:size(drawmapping_data,1)
            if all(drawmapping_data{mmhh,1} == this_monkey_hemi)
                xyoffsets = drawmapping_data{mmhh,2};      %  Changed from AP to refer to anterior commissure
                this_mapping_data = drawmapping_data{mmhh,3};
                
                for ee = 1:length(this_mapping_data)
                    %                         if all(this_mapping_data{ee}{2} == this_pos_raw(1:2))
                    
                    if all(this_mapping_data{ee}{1} == this_session) && all(this_mapping_data{ee}{2} == this_pos_raw(1:2))
                        found = found + 1;
                        
                        % Anterior-posterior
                        AP =  (this_pos_raw(1) - xyoffsets(1)) * 0.8; % Refer to AC0
                        %                             AP = - (this_pos_raw(1) - xyoffsets(1)) * 0.8; % in mm, P - --> AP0 --> + A
                        % Ventral-dorsal
                        %                             VD = (this_pos_raw(2) - xyoffsets(2)) * 0.8 - AP * tan(30/180*pi); % LIPs are inclined by about 30 degree
                        %                             % Commented by ZZ 20201229 for nearly no
                        %                             inclined. I'll adapt this value if it's
                        %                             indeed inclined in AP axis.
                        which_GM_is_CD = find((this_mapping_data{ee}{4}(:,1) == CD));
                        which_CD_this_cell_in =  find((this_mapping_data{ee}{4}(which_GM_is_CD,2) <= ceil(this_pos_raw(3)/100)) & ...
                            (this_mapping_data{ee}{4}(which_GM_is_CD,3) >= fix(this_pos_raw(3)/100)));
                        
                        if which_CD_this_cell_in == 1 % This cell is in the first (upper) layer of LIP
                            this_surface = this_mapping_data{ee}{4}(which_GM_is_CD(which_CD_this_cell_in),2); % Surface is the start of this GM
                            depth = this_pos_raw(3) - 100 * this_surface; % in um
                        elseif which_CD_this_cell_in == 2 % This cell is in the second (lower) layer of LIP, the depth should be inversed.
                            % this_surface = this_mapping_data{ee}{4}(which_GM_is_LIP(which_LIP_this_cell_in),3); % Surface is the end of this GM
                            % depth = - (this_pos_raw(3) - 100 * this_surface); % in um
                            depth = nan; % Because I'm not sure whether the end of penetration has reached the (bottom) surface of the second LIP
                        else
                            disp('Check CD layers!!!'); beep; keyboard
                        end
                        
                        %                             group_result(i).position = [this_monkey_hemi AP VD depth which_CD_this_cell_in]; % [monkey hemi AP Ventral-Dorsal depth which_LIP_this_cell_in]
                        group_result(i).position = [this_monkey_hemi AP depth which_CD_this_cell_in]; % [monkey hemi AP Ventral-Dorsal depth which_LIP_this_cell_in]
                        
                    end
                end
            end
        end
        
        if found~=1 % Should be exactly 1
            disp('Cell position not found!!');
            beep; keyboard;
        end
    end
    
    %         group_position = reshape([group_result.position],6,[])';
    group_position = reshape([group_result.position],5,[])';
    % changed by ZZ 20201229 for deleting VD iterm
    
    end


%% Final Preparation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===========  Common parameters   ============
% Overall time marker (Stim on, Stim off, Sac on)
% tt = reshape(cell2mat([group_result.time_marker]),6,[])';
tt = reshape(cell2mat([group_result.time_marker]),4,[])';  % I only have 4 time markers

time_markers{1} = [mean(tt(:,1:2)); std(tt(:,1:2))];
time_markers{2} = [mean(tt(:,3)); std(tt(:,3))];
time_markers{3} = [mean(tt(:,4)); std(tt(:,4))];

p_critical = 0.05;

% =========== Data for common use ============
function_handles = [];

% For PCA_A and hotgram
j_PCA_A = 1;
PCA_A = []; PCA_A_PC = []; sort_time_interval1 = []; sort_time_interval2 = []; sort_time_interval3 = [];
A_memSac = []; A_choicediv = []; A_moddiv = []; A_CP = [];
enlarge_factor = 30; % Enlarge memsac DDIs

% For PCA_B
j_PCA_B = 1;
select_for_PCA_B = select_bottom_line;
% PCA_B_time_range = min(rate_ts{j_PCA_B})+100 <= rate_ts{j_PCA_B} & rate_ts{j_PCA_B} <= time_markers{j_PCA_B}(1,3);  % Before saccade
PCA_B_time_range = min(rate_ts{j_PCA_B}) <= rate_ts{j_PCA_B} & rate_ts{j_PCA_B} <= max(rate_ts{j_PCA_B});  % Before saccade
% changed by ZZ 20201229

% PCA_B_time_range = min(rate_ts{j_PCA_B})+100 <= rate_ts{j_PCA_B} & rate_ts{j_PCA_B} <= time_markers{j_PCA_B}(1,2);  % Before stim
PCA_B_times = rate_ts{j_PCA_B}(PCA_B_time_range);
denoised_dim = 8;
PCA_B = []; weights_PCA_B_PC = []; PCA_B_projPC = []; PCA_B_explained = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ================ Miscellaneous ===================

% == Figure default
set(0, 'defaultFigureRenderer', 'painters')
set(0, 'defaultFigureRendererMode', 'Manual')
set(0,'defaultAxesColorOrder',[41 89 204; 248 28 83; 14 153 46]/255);
condition_colors = [41 89 204; 248 28 83; 14 153 46]/255;
heading_colors = [0 0 0; 0.3137 0.1961 0.1249; 0.6274 0.3921 0.2497; 0.9363 0.5851 0.3726; 1.0000 0.7812 0.4975];

modality_diff_colors = condition_colors;    modality_diff_colors(3,:) = [0 0 0];

transparent = get(findall(gcbf,'tag','transparent'),'value');
figN = 2499;

marker_for_time_markers{1} = {'-','-','--'};
marker_for_time_markers{2} = {'--','--','-'};


%% ====================================== Function Handles =============================================%

% Today I reconstruct the group analysis codes in a more logical and effcient way:
% 1. I divided each figure into small nested functions in Group_HD
% 2. Group_HD returns a cell called "function_handles" back to GROUP_GUI (with all the data available for each nested function)
% 3. Now I can choose which figure(s) to plot or which function to debug in GROUP_GUI
% 4. I feel so happy.
% @HH20150425

function_handles = {
    'Mean Rate Metrics', {
    'Correct only, all choices (conventional)', @f1p1;
    '   All single cell traces',@f1p1p6;
    '      Example cells', @f1p1p4;
    '      Divergence time', @f1p1p6p1;
    '   Different weighting methods',@f1p1p5;
    'Different headings' , @f1p2;
    '   Comb/max ratio distribution through time', @f1p2p5;
    '   Tuning curve analysis',@f1p2p7;
    '   Partial correlation (grand)',@f1p2p8;
    'Correct / Wrong Trials', @f1p3;
    };
    
    'ROC Metrics',{
    'CP, CDiv, and MDiv',@f2p1;
    'Multisensory Enhancement of CDiv',@f2p2;
    'Easy and Difficult',@f2p3;
    'Correct and Error', @f2p4;
    };
    
    'Correlations', {
    'Mem-sac vs. pre/post CDiv/CP, CDiv vs. CP', @f3p1;
    'Mem-sac vs. abs(CPref)', @f3p1p2;
    'Choice Preference vs. Modality Preference', @f3p2;
    'Choice Preference between modalities', @f3p2p2;
    'Choice Preference: pre and post', @f3p2p3;
    'Choice Preference: Correct and Error', @f3p2p4;
    'ChoiceDiv: Correct and Error', @f3p2p5;
    'Choice Preference: Difficult and Easy', @f3p2p6;
    'ChoiceDiv: Difficult and Easy', @f3p2p7;
    'Psychophysics vs. CD', @f3p3;
    'Cell position Distribution',@f3p4;
    };
    
    'PCA_choice & modality (Eigen-feature)',{
    'Hot-gram', @f4p1;
    'Cluster and Trajectory', @f4p2;
    };
    
    'PCA_choice & modality (Eigen-neuron)',{
    'Weights and correlations', @f5p1;
    '1-D Trajectory',@f5p2;
    '3-D Trajectory',@f5p3;
    };
    
    'PCA_heading & choice (Eigen-neuron)',{
    '2-D Trajectory', @f51p1;
    };
    
    'dPCA & TDR', {
    'demixed PCA', @f5p4;
    'TDR_ZZ', @f5p5; 
    }
    
    'Linear SVM decoder (choice and modality)',{
    'Training SVM', @f6p0;
    'Weights', @f6p1;
    'Performance (overall)', @f6p2;
    'SVM weighted sum',@f6p3;
    };
    
    'Lasso Decoder', {
    'Data Preparing', @f11p0;
    'Model Training', @f11p1;
    'Visualizing Model', @f11p2;
    'Cross-modal Testing', @f11p3;
    'Cross-modal Accuracy', @f11p4;
    'Within-modal Testing', @f11p5;
    'Within-modal Accuracy', @f11p6;
    };
    
        'Fisher information of heading',{
    'Simple Fisher like Gu 2010: Sum(slope/mean)', @f6p5p1
    '   + Correlation between modalities',@f6p5p1p1
    '   + Correlation with CPref',@f6p5p1p2
    'Dora''s partial corr and linear regression k', @f6p5p2
    'Training SVM decoders', @f6p5p9;
    }
    
    'Targeted Dimensionality Reduction',{
    'PCA + SVM: weights', @f8p1;
    'PCA + SVM: all correct + all angles',@f8p2;
    }
    
    'Others',{
    'Cell Counter',@f9p1;
    'Target first vs. Target last',@f9p2;
    '-----------------------------------', @cell_selection;
    'Export Associated Memsac Files',@f9p9;
    };
    
    'NoShow',{@cell_selection};
    
    };

%% ====================================== Function Definitions =============================================%

    function f1p1(debug)      % Rate 1. Correct only, all choices
        if debug  ; dbstack;   keyboard;      end
        
        %% ------- Averaged norm PSTH --------
        set(figure(999),'name','Average PSTH (Correct only, all choices)','pos',[27 63 1449 892]); clf
        h_subplot = tight_subplot(2,3,[0.11 0.05],[0.05 0.1],0.07);
        
        methods_of_select = {
            select_bottom_line, 'All cells';
            select_tcells, 'Typical cells';
            select_no_tcells, 'Non-typical cells'};
        
        for j= 1:3
            for ms = 1:3
                % --- Difference (Pref - Null) ---
                PSTH_all_Norm_PrefminusNull{j}{ms} = PSTH_all_Norm{j}(methods_of_select{ms,1},:,1:2:end)...
                    - PSTH_all_Norm{j}(methods_of_select{ms,1},:,2:2:end);
            end
        end
        
        for ms = 1:3
            
            h = SeriesComparison({PSTH_all_Norm{1}(methods_of_select{ms,1},:,:), PSTH_all_Norm{2}(methods_of_select{ms,1},:,:), PSTH_all_Norm{3}(methods_of_select{ms,1},:,:)},...
                {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
                'Colors',{condition_colors(1,:),condition_colors(1,:),condition_colors(2,:),condition_colors(2,:),condition_colors(3,:),condition_colors(3,:)},'LineStyles',{'-','--','-','--','-','--'},...
                'ErrorBar',6,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot((ms-1)*2+1),...
                'CompareIndex',[1,3,5;2,4,6],...
                'CompareColor',[mat2cell(condition_colors,ones(3,1))],...
                'Transparent',transparent, 'PCritical',0.01);
            
            xlabel('Time (ms)');
            legend off;
            title([methods_of_select{ms,2} ', n = ' num2str(sum(methods_of_select{ms,1}))]);
            axis tight;
            
            %             xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);
%             ylim([0.05 0.6]);
            
            % Gaussian vel
            plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            
            %             legend off;
            
            SeriesComparison({PSTH_all_Norm_PrefminusNull{1}{ms}, PSTH_all_Norm_PrefminusNull{2}{ms},PSTH_all_Norm_PrefminusNull{3}{ms}},...
                {rate_ts{1} rate_ts{2}  rate_ts{3} time_markers},...
                'Colors',mat2cell(condition_colors,ones(3,1)),'LineStyles',{'-'},...
                'ErrorBar',6,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot((ms-1)*2+2),...
                'CompareIndex',[1:3;1:3],...   % Change from [1:3,1,2; 1:3, 3,3] to this ZZ 20210517
                'CompareColor',[mat2cell(condition_colors,ones(3,1));condition_colors(1,:);condition_colors(2,:);condition_colors(3,:)],...
                'Transparent',transparent, 'PCritical',0.01);
            
            
            xlabel('Time (ms)');
            legend off;
            
            axis tight;
            %             xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);
%             ylim([-0.05 0.35]);
            
            % Gaussian vel
            plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/6,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            
            legend off;
            
        end
        %         end
        SetFigure(15);
        
        %% ------- Averaged raw PSTH --------
        set(figure(1999),'name','Average PSTH (Correct only, all choices)','pos',[27 63 1449 892]); clf
        h_subplot = tight_subplot(2,3,[0.11 0.05],[0.05 0.1],0.07);
        
        methods_of_select = {
            select_bottom_line, 'All cells';
            select_tcells, 'Typical cells';
            select_no_tcells, 'Non-typical cells'};
        
        for j = 1:3
            for ms = 1:3
                PSTH_all_raw_PrefminusNull_tmp{j}{ms} = PSTH_all_raw{j}(methods_of_select{ms,1},:,1:2:end)...
                    - PSTH_all_raw{j}(methods_of_select{ms,1},:,2:2:end);
                
            end
        end
        for ms = 1:3
            
            h = SeriesComparison({PSTH_all_raw{1}(methods_of_select{ms,1},:,:), PSTH_all_raw{2}(methods_of_select{ms,1},:,:),PSTH_all_raw{3}(methods_of_select{ms,1},:,:)},...
                {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
                'Colors',{condition_colors(1,:),condition_colors(1,:),condition_colors(2,:),condition_colors(2,:),condition_colors(3,:),condition_colors(3,:)},'LineStyles',{'-','--','-','--','-','--'},...
                'ErrorBar',6,'Xlabel',[],'Ylabel','Raw firing','axes',h_subplot((ms-1)*2+1),...
                'CompareIndex',[1,3,5; 2,4,6],...
                'CompareColor',[mat2cell(condition_colors,ones(3,1))],...
                'Transparent',transparent, 'PCritical',0.01);
            
            %             if ms < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
            
            xlabel('Time (ms)');
            legend off;
            title([methods_of_select{ms,2} ', n = ' num2str(sum(methods_of_select{ms,1}))]);
            axis tight;
            
            %             xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);
%             ylim([0 15]);
            
            % Gaussian vel
            plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            
            
            SeriesComparison({PSTH_all_raw_PrefminusNull_tmp{1}{ms}, PSTH_all_raw_PrefminusNull_tmp{2}{ms}, PSTH_all_raw_PrefminusNull_tmp{3}{ms}},...
                {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
                'Colors',mat2cell(condition_colors,ones(3,1)),'LineStyles',{'-'},...
                'ErrorBar',6,'Xlabel',[],'Ylabel','Raw firing','axes',h_subplot((ms-1)*2+2),...
                'CompareIndex',[1:3;1:3],...
                'CompareColor',[mat2cell(condition_colors,ones(3,1));condition_colors(1,:);condition_colors(2,:);condition_colors(3,:)],...
                'Transparent',transparent, 'PCritical',0.01);
            
            
            xlabel('Time (ms)');
            legend off;
            axis tight;
            %             xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);
%             ylim([-1 10]);
            
            % Gaussian vel
            plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/6,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            
        end
        SetFigure(15);
        
    end
    function f1p1p5(debug)
        if debug  ; dbstack;   keyboard;      end
        
        % ------- Weighted by different methods ------
        % Straight mean
        weights_straight_mean = ones(sum(select_bottom_line),1);
        Weighted_sum_PSTH( weights_straight_mean,{'Weighted by straight mean'},select_bottom_line);
        title('Weighted by same weights (straight mean)');
        Weights_Property_Correlation( weights_straight_mean /sum(weights_straight_mean) ,...
            {'Weighted by straight mean','Weighted by straight mean'},select_bottom_line);
        
        % Norm by dynamic range
        weights_norm = weights_normalized_PSTH(select_bottom_line)';
        Weighted_sum_PSTH( weights_norm,{'Weighted by 1/dynamic range'},select_bottom_line);
        title('Weighted by 1/dynamic range');
        Weights_Property_Correlation(weights_norm/sum(weights_norm),...
            {'Weighted by 1/dynamic range','Weighted by 1/dynamic range'},select_bottom_line);
        
        % Norm by dynamic range with cell selection
        weights_norm_tcells = weights_norm;
        weights_norm_tcells(~select_tcells(select_bottom_line)) = 0;
        Weighted_sum_PSTH( weights_norm_tcells,{'Weighted by 1/dynamic range with cell selection'},select_bottom_line);
        title('Weighted by 1/dynamic range with cell selection');
        Weights_Property_Correlation(weights_norm_tcells/sum(weights_norm_tcells),...
            {'Weighted by 1/dynamic range with cell selection','Weighted by 1/dynamic range with cell selection'},select_bottom_line);
        
    end

%% =================================================================================================================
% For f1p1p6 & f1p1p7. @HH20160906
% All angles, Pref - null
PSTH_all_raw_PrefminusNull{1} = PSTH_all_raw{1}(select_bottom_line,:,1:2:end)...
    - PSTH_all_raw{1}(select_bottom_line,:,2:2:end);
PSTH_all_raw_PrefminusNull{2} = PSTH_all_raw{2}(select_bottom_line,:,1:2:end)...
    - PSTH_all_raw{2}(select_bottom_line,:,2:2:end);
PSTH_all_raw_PrefminusNull{3} = PSTH_all_raw{3}(select_bottom_line,:,1:2:end)...
    - PSTH_all_raw{3}(select_bottom_line,:,2:2:end);
[single_cell_plot_order_index,single_cell_plot_order] = sort(mean(abs(Choice_pref_all(:,select_bottom_line(select_cpref_mpref),1))),'descend'); % Stim-on to stim-off

    function f1p1p6(debug)
        if debug  ; dbstack;   keyboard;      end
        %%
        set(figure(700),'name','Single cell (Correct only, all choices)','pos',[27 63 1449 892]); clf
        [~,h_subplot] = tight_subplot(4,5,[0.03 0.02]);
        counter = 0;
        ps = nan(3,length(rate_ts{1}),length(single_cell_plot_order));
        in_larger_than_out = ps;
        
        for nn = 1:length(single_cell_plot_order)
            counter = counter + 1;
            
            if counter > 20
                counter = 1;
                set(figure(700+fix(nn/20)),'name','Single cell (Correct only, all choices)','pos',[27 63 1449 892]); clf
                [~,h_subplot] = tight_subplot(4,5,[0.03 0.02]);
            end
            
            cellNo_in_find_bottom_line = single_cell_plot_order(nn);
            cellNo_origin = find(cumsum(select_bottom_line)==cellNo_in_find_bottom_line,1);
            
            for j = 1:3
                
                ys_this{j} = group_result(cellNo_origin).mat_raw_PSTH.PSTH{j,1,1}.ys';
                sem_this{j} = group_result(cellNo_origin).mat_raw_PSTH.PSTH{j,1,1}.sem';
                ps_this{j} = group_result(cellNo_origin).mat_raw_PSTH.PSTH{j,1,1}.ps';
            end
            
            SeriesComparison({shiftdim(ys_this{1},-1) shiftdim(ys_this{2},-1) shiftdim(ys_this{3},-1)},...
                {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
                'OverrideError',{sem_this{1}, sem_this{2} sem_this{3}},...
                'OverridePs',{ps_this{1}, ps_this{2} ps_this{3}},'ErrorBar',6,...
                'CompareIndex',[1 3 5;2 4 6],'CompareColor',{condition_colors(1,:),condition_colors(2,:),[0 0.8 0.4]},...
                'Colors',{condition_colors(1,:),condition_colors(1,:),condition_colors(2,:),condition_colors(2,:),[0 0.8 0.4],[0 0.8 0.4]},'LineStyles',{'-','--'},'axes',h_subplot(counter),...
                'Transparent',transparent, 'PCritical',0.01);
            
            title(h_subplot(counter), sprintf(' %g (ori. #%g, ps: %g, %g, %g; %g)',nn,cellNo_origin,Choice_pref_p_value_all(:,cellNo_in_find_bottom_line,1)',single_cell_plot_order_index(nn)));  % 3, stim-on to stim-off
            legend off;
            axis(h_subplot(counter),'tight');
            ylabel(h_subplot(counter), '');
            
            % Cache for calcuating the first significant times. HH20161108
            ps(:,:,nn) = ps_this{1}';
            in_larger_than_out(:,:,nn) = (ys_this{1}(:,1:2:end) > ys_this{1}(:,2:2:end))';
            
            % Annotate typical cell
            if select_tcells(cellNo_origin)
                set(gca,'color',hsv2rgb([1 0.1 1]));
            end
            
            % Plot mem-sac trace
            set(gca,'ButtonDownFcn',{@Plot_HD,cellNo_origin});
            
            if counter ~= 20
                set(h_subplot(counter),'xticklabel','');
            end
            
            plot(h_subplot(counter),xlim,[0 0],'k--');
            
        end
        %}
    end

firstDivergenceTime = [];

    function f1p1p6p1(debug)  % Divergence time
        if debug  ; dbstack;   keyboard;      end
        
        %% === Divergence time ===
        %%%%%%%%%%%%%%%%%%%%%%%
        min_bins_for_1st_sign_time = 30; % 300 ms
        p_critical_for_1st_sign_time = 0.05;    % Paired t-test for two line during each time window
        t_valid_ind = rate_ts{1}>= 200 & rate_ts{1}<=1500; % Only sensory phase
        %%%%%%%%%%%%%%%%%%%%%%%
        
        firstDivergenceTime = nan(sum(select_bottom_line),3);
        t_valid = rate_ts{1}(t_valid_ind);
        
        monkeys = xls_num{1}(:,header.Monkey);  % HH20160613
        monkeys = monkeys(select_bottom_line); 
        %         monkey1 = monkeys == 15; monkey1 = monkey1(select_bottom_line)';
        %         monkey2 = monkeys == 13; monkey2 = monkey2(select_bottom_line)';
        
        cell_ind = find(select_bottom_line); 

        for nn = 1:sum(select_bottom_line)
            ys_this = group_result(cell_ind(nn)).mat_raw_PSTH.PSTH{1,1,1}.ys'; % j = 1
            ps_this = group_result(cell_ind(nn)).mat_raw_PSTH.PSTH{1,1,1}.ps;  % j = 1
            in_larger_than_out = (ys_this(:,1:2:end) > ys_this(:,2:2:end))';
            
            ps_valid = double((ps_this(:,t_valid_ind) < p_critical_for_1st_sign_time) & in_larger_than_out(:,t_valid_ind));
            
            conv_tmp = conv2(ps_valid,ones(1,min_bins_for_1st_sign_time));
            for k = 1:3
                % This is really brilliant :)
                first_tmp = find(conv_tmp(k,:) == min_bins_for_1st_sign_time,1) - min_bins_for_1st_sign_time + 1;
                
                if ~isempty(first_tmp)
                    firstDivergenceTime(nn,k) = t_valid(first_tmp);
                end
            end
        end
        
        %         first_sign_result = BarComparison(firstDivergenceTime,'figN',556,'Colors',mat2cell(colors,ones(stim_type_num,1)));
        %         first_sign_ps =  [first_sign_result.ps_ttest(2,3),first_sign_result.ps_ttest(2,4),first_sign_result.ps_ttest(3,4)];
        %         title(sprintf('%g, ',sum(~isnan(firstDivergenceTime)),first_sign_ps));
        %         ylim([0 1500]); SetFigure(15);
        %         view([90 90]);
        set(figure(figN), 'pos',[50 50 1000 800], 'Name','Divergence Time'); clf; hold on; figN=figN+1;
        
        for m = 1:length(monkey_included_for_analysis)
            %             for m = monkey_included_for_analysis
            %             m = find(monkey_included_for_analysis == monkey_included_for_loading);
            mean_firstDiverTime(m,:) = nanmean(firstDivergenceTime(monkeys==monkey_included_for_analysis(m),:),1);
            stds_firstDiverTime(m,:) = nanstd(firstDivergenceTime(monkeys==monkey_included_for_analysis(m),:),1);
            sems_firstDiverTime(m,:) = stds_firstDiverTime(m,:)./sqrt(sum(~isnan(firstDivergenceTime(monkeys==monkey_included_for_analysis(m),:)),1));
            
            % Drop the Bar
            % ZZ @ 20220308
            plot((1:3) - 0.05+0.05*m, firstDivergenceTime(monkeys==monkey_included_for_analysis(m),:),monkey_line{m},'Color',[0.8 0.8 0.8], 'LineWidth',1.5); hold on;
            
            for c = 1:stim_type_num
                scatter((c-0.1+0.1*m)*(~isnan(firstDivergenceTime(monkeys==monkey_included_for_analysis(m),c))),...
                    firstDivergenceTime(monkeys==monkey_included_for_analysis(m),c),...
                    150, monkey_marker{m},...
                    'MarkerEdgeColor','k', 'MarkerFaceColor',condition_colors(c,:),...
                    'MarkerEdgeAlpha', 0.5, 'MarkerFaceAlpha',0.5);
                
                yneg = [sems_firstDiverTime(m,c)];
                ypos = [sems_firstDiverTime(m,c)];
                
                errorbar(c-0.25+0.2*m, mean_firstDiverTime(m,c),yneg, ypos,monkey_marker{m},'Color',condition_colors(c,:),...
                    'MarkerEdgeColor',condition_colors(c,:),'MarkerFaceColor',condition_colors(c,:), 'MarkerSize',12,'CapSize',18,...
                    'LineWidth', 2, 'CapSize',20);
                
            end
            
        end
        
        % Paired t-test
        for c = 1:stim_type_num
            for c2 = c+1:size(firstDivergenceTime,2)
                
                if length(monkey_included_for_analysis)>1
                    [~,p] = ttest(firstDivergenceTime(:,c), firstDivergenceTime(:,c2));  % test as a whole
                    ps_ttest(c,c2) = p;
                    
                    text((c+c2)/2-0.2, 1400, {['n= ' num2str(sum(all(~isnan(firstDivergenceTime(:,[c c2])),2)))],...
                        ['p = ', num2str(p)]});
                    
                else
                    [~,p] = ttest(firstDivergenceTime(monkeys==monkey_included_for_analysis(m),c), firstDivergenceTime(monkeys==monkey_included_for_analysis(m),c2));
                    %                 [~,p] = ttest2(firstDivergenceTime(:,c), firstDivergenceTime(:,c2));
                    ps_ttest(c,c2) = p;
                    
                    text((c+c2)/2-0.2, 1400, {['n= ' num2str(sum(all(~isnan(firstDivergenceTime(monkeys==monkey_included_for_analysis(m),[c c2])),2)))],...
                        ['p = ', num2str(p)]});
                    
                end
                
            end
        end
        
        ylim([0 1500]); xlim([0.8 3.2]);
        xticks(1:3);
        
        text(1.5, 10, num2str(mean_firstDiverTime))
        text(2.5, 10, num2str(sems_firstDiverTime))
        
        SetFigure();
        view([90 90]);
    end


selected_t_vest;
selected_t_vis;
selected_t_comb;

    function f1p2(debug)      % Rate 2. Different headings
        %% Different headings
        
        if debug
            dbstack;
            keyboard;
        end
        
        %         methods_of_select = {
        %             % select_bottom_line, 'All cells';
        %             select_tcells, 'Typical cells';
        %             };
        
        methods_of_select = {
            % select_bottom_line, 'All cells';
            selected_t_vest, 'Vest. cell';
            selected_t_vis, 'Vest. cell';
            selected_t_comb, 'Vest. cell';
            select_tcells, 'Typical cells';
            };
        
        colors_angles = reshape(repmat(heading_colors, 1, 2)',3, [])';
        colors_angles = mat2cell(colors_angles, ones(10,1));
        
        for ms = 1 % :size(methods_of_select,1)
            set(figure(999-ms),'name',['Average PSTH (correct only, all headings), ' methods_of_select{ms,2}],'pos',[27 57 919 898]); clf
            h_subplot = tight_subplot(3,2,[0.05 0.1],[0.05 0.15],[0.12 0.03]);
            %             linkaxes(h_subplot(1:3),'xy')
            %             linkaxes(h_subplot(4:6),'xy')
            
            for k = 1:3
                
                %                 colors_angles = colormap(gray);
                %                 colors_angles = ones(size(colors_angles,1),3) - colors_angles .* repmat([1 1 1]-colors(k,:),size(colors_angles,1),1);
                %                 colors_angles = colors_angles(round(linspace(20,length(colors_angles),5)),:);
                %                 colors_angles = reshape(repmat(colors_angles,1,2)',3,[])';
                %                 colors_angles = mat2cell(colors_angles,ones(10,1));
                
                % --- Ramping with different angles ---
                for j = 1:3
                    %                      yyy{j} = PSTH_correct_angles_Norm{j}(methods_of_select{ms,1},:,:,k);
                    yyy{j} = PSTH_correct_angles_raw{j}(methods_of_select{k,1},:,:,k);
                    ttt{j} = rate_ts{j};
                    
                    yyy_diff{k}{j} =  yyy{j}(:,:,1:2:end) - yyy{j}(:,:,2:2:end);
                end
                
                
                h = SeriesComparison(yyy,{ttt{1} ttt{2} ttt{3} time_markers},...
                    'Colors',colors_angles,'LineStyles',{'-','--'},...
                    'ErrorBar',0,'Xlabel',[],'Ylabel','Raw firing','axes',h_subplot(k));
                
                if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                
                %                 legend off;
                
%                 if k == 1
                    title([methods_of_select{k,2} ', n = ' num2str(sum(methods_of_select{k,1}))]);
%                 end
                
                xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);  % ylim([0.1 max(ylim)]);
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                
                legend off;
                
                
                % ----  Difference ---
                
                h = SeriesComparison(yyy_diff{k},{ttt{1} ttt{2} ttt{3} time_markers},...
                    'Colors',colors_angles(1:2:end),'LineStyles',{'-'},...
                    'ErrorBar',0,'Xlabel',[],'Ylabel','Diff','axes',h_subplot(k+(2-1)*3));
                
                if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                
                %                 legend off;
                
                xlim([rate_ts{1}(10) rate_ts{1}(end-10)]); ylim(ylim); plot(xlim,[0 0],'k--');
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                
                legend off;
                
            end
            
            %             axis(h_subplot(end),'tight');
            SetFigure(15); drawnow;
            
            
            % --- Correct and Wrong trials of different angles --- HH20180711 for GuYong's Hangzhou ppt
            
            set(figure(221816-ms),'name',['Average PSTH (correct + wrong, all headings), ' methods_of_select{ms,2}]); clf
            h_subplot = tight_subplot(3,3,[0.05 0.1],[0.05 0.15],[0.12 0.03]);
            set(gcf,'uni','norm','pos',[0.009       0.059       0.732       0.853]);
            %             linkaxes(h_subplot,'xy')
            
            for k = 1:3
                
                %                 colors_angles = colormap(gray);
                %                 colors_angles = ones(size(colors_angles,1),3) - colors_angles .* repmat([1 1 1]-colors(k,:),size(colors_angles,1),1);
                %                 colors_angles = colors_angles(round(linspace(20,length(colors_angles),5)),:);
                %                 colors_angles = reshape(repmat(colors_angles,1,2)',3,[])';
                %                 colors_angles = mat2cell(colors_angles,ones(10,1));
                
                % --- Ramping with different angles ---
                for j = 1:3
                    yyy{j} = PSTH_correct_angles_Norm{j}(methods_of_select{k,1},:,:,k);
                    %                     yyy{j} = PSTH_correct_angles_raw{j}(methods_of_select{ms,1},:,:,k);
                    ttt{j} = rate_ts{j};
                    
                    yyy_diff{k}{j} =  yyy{j}(:,:,1:2:end) - yyy{j}(:,:,2:2:end);
                end
                
                h = SeriesComparison(yyy,{ttt{1} ttt{2} ttt{3} time_markers},...
                    'Colors',colors_angles,'LineStyles',{'-','--'},...
                    'ErrorBar',0,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(k));
                
                if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                
                legend off;
                
%                 if k == 1
                    title([methods_of_select{k,2} ', Correct only, n = ' num2str(sum(methods_of_select{k,1}))]);
%                 end
                
                xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);  % ylim([0.1 max(ylim)]);
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                
                legend off;
                
                
                % --- Ramping with different angles, Correct + Wrong ---
                % The order has been aligned, from non-PREF to PREF heading
                % But not the raw headings
%                 colors_angles = mat2cell([flipud(heading_colors(2:end,:)); heading_colors], ones(9,1));
                
                
                %                 colors_angles = colormap(gray);
                %                 colors_angles = ones(size(colors_angles,1),3) - colors_angles .* repmat([1 1 1]-colors(k,:),size(colors_angles,1),1);
                %                 colors_angles = colors_angles(round(linspace(20,length(colors_angles),5)),:);
                %                 colors_angles = [flipud(colors_angles); colors_angles(2:end,:)];
                
                for j = 1:3
                    yyy{j} = PSTH_correctNwrong_angles_Norm{j}(methods_of_select{k,1},:,:,k);
                    %                     yyy{j} = PSTH_correctNwrong_angles_raw{j}(methods_of_select{ms,1},:,:,k);
                end
                
                h = SeriesComparison(yyy,{rate_ts{1} rate_ts{2} rate_ts{2} time_markers},...
                    'Colors',[colors_angles(9:-2:1); colors_angles(3:2:end)],'LineStyles',{'--' ,'--','--' ,'--' ,'-', '-', '-', '-', '-'},...
                    'ErrorBar',0,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(k+(2-1)*3));
                
                if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                
                legend off;
                
%                 if k == 1
                    title([methods_of_select{k,2} ', Correct + Wrong, n = ' num2str(sum(methods_of_select{k,1}))]);
%                 end
                
                xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);  %ylim([0.1 0.8]);
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                
                %{
                % --- Ramping with different angles, Wrong only ---
                % Has been changed to [small pref, small null, ..., large pref, large null]
                
                for j = 1:3
                    yyy{j} = PSTH_wrong_angles_Norm{j}(methods_of_select{ms,1},:,:,k);
                    %                     yyy{j} = PSTH_wrong_angles_raw{j}(methods_of_select{ms,1},:,:,k);
                end
                
                % Note that for wrong only, here I align to the "stimulus" not the "choice"
                h = SeriesComparison(yyy,{rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
                    'Colors',colors_angles([4 4 3 3 2 2 1 1],:),'LineStyles',{'--','-'},...      %  dashed line indicate choose Non-PREF direction
                    'ErrorBar',0,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(k+(3-1)*3));
                
                if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                
                legend off;
                
                if k == 1
                    title([methods_of_select{ms,2} ', Wrong only, n = ' num2str(sum(methods_of_select{ms,1}))]);
                end
                
                xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);  %ylim([0.1 0.8]);
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                
                
                legend off;
                
                %}
            end
            
            SetFigure(15); drawnow;
            %             axis(h_subplot(3),'tight');
            
%             ms = 4; 
%             % --- Multisensory enhancement of different angles --- HH20160126
%             set(figure(995-ms),'name',['Enhancement for different angles (correct only, all headings), ' methods_of_select{ms,2}],'pos',[218 35 1674 928]); clf
%             h_subplot = tight_subplot(2,3,[0.1 0.05],[0.05 0.07],[0.1 0.05]);
%             h_subplot =  reshape(reshape(h_subplot,2,3)',[],1); delete(h_subplot(end));
%             unique_abs_heading = unique(abs(group_result(representative_cell).mat_raw_PSTH.heading_per_trial));
%             %             linkaxes(h_subplot,'xy')
%             
%             for aa = 1:size(yyy_diff{1}{1},3) % Different angles
%                 
%                 % === Plotting ===
%                 yyy_diff_this_angle_all_stim_type = [];
%                 for j = 1:3
%                     for k = 1:3
%                         yyy_diff_this_angle_all_stim_type{j}(:,:,k) = yyy_diff{k}{j}(:,:,aa); % Reorganize data
%                     end
%                 end
%                 
%                 h = SeriesComparison(yyy_diff_this_angle_all_stim_type, {ttt{1} ttt{2} ttt{2} time_markers},...
%                     'Colors',mat2cell(condition_colors,ones(3,1)),'LineStyles',{'-'},...
%                     'ErrorBar',2,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(aa),'Transparent',0);
%                 xlim([rate_ts{1}(10) rate_ts{1}(end-10)]);
%                 %                 ylim([-.1 .5]);
%                 plot(xlim,[0 0],'k--');
%                 
%                 
%                 title(['|heading| = ',num2str(unique_abs_heading(aa))]);
%                 if k == 1
%                     title([methods_of_select{ms,2} ', n = ' num2str(sum(methods_of_select{ms,1}))]);
%                 end
%                 
%                 if aa ~= 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
%                 
%                 % Gaussian vel
%                 plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
%                 legend off;     drawnow
%                 
%                 % === Linear fitting combined trace. HH20160510 ===
%                 ts = rate_ts{1};
%                 t_select = rate_ts{1}> 0 & rate_ts{1} <=1500;
%                 
%                 ramping1(:,aa) = h.means{1}(1,t_select);
%                 ramping2(:,aa) = h.means{1}(2,t_select);
%                 ramping3(:,aa) = h.means{1}(3,t_select);
%                 
%                 w = fminsearch(@(w) sum((w(1)*ramping1(:,aa) + w(2)*ramping2(:,aa) - ramping3(:,aa)).^2), [.5 .5]);
%                 
%                 plot(ts(t_select),ramping1(:,aa)*w(1)+ramping2(:,aa)*w(2),'k','linew',2);
%                 %                 plot(ts(t_select),ramping1(:,aa)*sqrt(2)/2+ramping2(:,aa)*sqrt(2)/2,'m--','linew',2);
%                 text(700,0,num2str(w));
%             end
%             %             axis(h_subplot(end-1),'tight');
%             
%             % == All angle use the same weight. HH20160926 ==
%             w_all = fminsearch(@(w) sum((w(1)*ramping1(:) + w(2)*ramping2(:) - ramping3(:)).^2), [.5 .5]);
%             for aa = 1:size(yyy_diff{1}{1},3) % Different angles
%                 axes(h_subplot(aa));
%                 plot(ts(t_select),ramping1(:,aa)*w_all(1)+ramping2(:,aa)*w_all(2),'c','linew',2);
%                 text(700,0.1,num2str(w_all),'color','c');
%             end
%             
%             SetFigure(15);
        end
        
    end

    function f1p1p4(debug) % Plot Example cells for Supplementary Figure 2
        if debug  ; dbstack;   keyboard;      end
        
        to_plot_ori = {
            % Four example cells in Figure 1
            %{
            [
            61 % m5c313r2u05
            126 % m10c133r3_4_7
            31 % m5c174r1
            120
            ]; % m10c121r3
            %}
            % Other exmples in SuppleFig 2
            %         %{
            [
            104
            106
            112
            121
            129
            131
            136
            
            
            ]
            %}
            [  % Bad cells
            171
            179
            184
            185
            199
            202
            205
            ]
            
            
            }; % Manually selected
        
        toPlotFigN = length(to_plot_ori);
        for ff = 1:toPlotFigN
            
            toPlotN = length(to_plot_ori{ff});
            
            figure(5837+ff); clf; % A large figure
            set(gcf,'uni','norm','pos',[0.001   0.04    0.626    0.78/4*toPlotN]);
            hs = tight_subplot(toPlotN,5,[0.06*4/toPlotN 0.06],[0.1 0.1]*4/toPlotN, [0.05 0.05]);
            
            hs = reshape(hs,[],5);
            
            for ee = 1:toPlotN
                Plot_HD([],[],to_plot_ori{ff}(ee),hs(ee,:));
            end
            
            % set(hs(1:end-1,:),'xtick','')
        end
    end


    function f1p2p5(debug)      % Rate 2.5. Different headings, Gu's idea, Enhance index (comb/max ratio) distribution through time
        if debug
            dbstack;
            keyboard;
        end
        
        methods_of_select = {
            select_tcells, 'Typical cells';
            };
        
        unique_abs_heading = unique(abs(group_result(representative_cell).mat_raw_PSTH.heading_per_trial));
        
        for ms = 1:1
            
            j = 1;
            ttt = rate_ts{j};
            
            % Organize data
            yyy_diff_angle_cellBtBheadingBk{j} =  PSTH_correct_angles_raw{j}(methods_of_select{ms,1},:,1:2:end,:)...
                - PSTH_correct_angles_raw{j}(methods_of_select{ms,1},:,2:2:end,:);
            yyy_diff_hard_easy_cellBtB2Bk{j} = PSTH_hard_easy_raw_cellBtB4Bk{j}(methods_of_select{ms,1},:,1:2:end,:)...
                - PSTH_hard_easy_raw_cellBtB4Bk{j}(methods_of_select{ms,1},:,2:2:end,:);
            
            % Sliding window
            sliding_width = 300; % ms, overlap allowed
            sliding_step = 200;
            ttt_1 = -100;
            sliding_n = ceil((ttt(end)-ttt_1-sliding_width)/sliding_step);
            
            % Hist range
            hist_range = (-1:0.1:4) + 0.05;
            plot_range = [-1.1 4.15];
            
            % -- Grouped by headings --
            for aa = 1:size(yyy_diff_angle_cellBtBheadingBk{j},3) % Different angles
                
                set(figure(700+aa),'name',['Enhance index, heading = ' num2str(unique_abs_heading(aa)) ', ' methods_of_select{ms,2}],'pos',[2 379 1674 581]); clf
                h_subplot = tight_subplot(1,sliding_n,[0,0.02],[0.2,0.1]);
                maxmaxy = 0;
                
                for ss = 1:sliding_n % Sliding window
                    
                    EI_t_range = (ss-1)*sliding_step + ttt_1 <= ttt & ttt <= (ss-1)*sliding_step + ttt_1 + sliding_width;
                    EI_ts = mean(ttt(EI_t_range));
                    
                    yyy_diff_this_angle_this_bin_cellBstimtype = ...
                        squeeze(mean(yyy_diff_angle_cellBtBheadingBk{j}(:,EI_t_range,aa,:),2));
                    
                    EI_this_bin = yyy_diff_this_angle_this_bin_cellBstimtype(:,3)./...,
                        max(yyy_diff_this_angle_this_bin_cellBstimtype(:,1),...
                        yyy_diff_this_angle_this_bin_cellBstimtype(:,2));
                    %                     EI_this_bin = yyy_diff_this_angle_this_bin_cell_by_stimtype(:,3)./...,
                    %                         (yyy_diff_this_angle_this_bin_cell_by_stimtype(:,2)+...
                    %                         0*yyy_diff_this_angle_this_bin_cell_by_stimtype(:,2));
                    
                    median_EI_this_bin = nanmedian(EI_this_bin);
                    mean_EI_this_bin = nanmean(EI_this_bin);
                    [~,p] = ttest(EI_this_bin,1);
                    
                    axes(h_subplot(ss));
                    hist(EI_this_bin,hist_range);    hold on; axis tight;
                    view(270,90); xlim(plot_range);
                    
                    set(gca,'xaxislocation','top');
                    
                    if ss<sliding_n
                        %set(gca,'xticklabel',[],'yticklabel',[]);
                        set(gca,'xticklabel',[]);
                    end
                    
                    ylabel(sprintf('%g ~ %g',ttt(find(EI_t_range,1)),ttt(find(EI_t_range,1,'last'))));
                    title(sprintf('med %.3g\nmean %.3g\np = %.2g',median_EI_this_bin,mean_EI_this_bin,p));
                    
                    maxmaxy = max(maxmaxy,max(ylim));
                    
                    plot([median_EI_this_bin median_EI_this_bin],[0 20],'r-','linew',2);
                    plot([mean_EI_this_bin mean_EI_this_bin],[0 20],'g-','linew',2);
                    
                    % Save EI
                    EI_all_angles_cellBtimeBheading(:,ss,aa) = EI_this_bin;
                end
                
                for ss = 1:sliding_n
                    %                     ylim(h_subplot(ss),[0 maxmaxy]);
                end
                
                SetFigure(15);
                
            end
            
            % -- Grouped by hard/easy --
            for dd = 1:2
                maxmaxy = 0;
                set(figure(750+dd),'name',['Enhance index, [hard, easy] = ' num2str(dd) ', ' methods_of_select{ms,2}],'pos',[2 379 1674 581]); clf
                h_subplot = tight_subplot(1,sliding_n,[0,0.02],[0.2,0.1]);
                
                for ss = 1:sliding_n % Sliding window
                    
                    EI_t_range = (ss-1)*sliding_step + ttt_1 <= ttt & ttt <= (ss-1)*sliding_step + ttt_1 + sliding_width;
                    
                    yyy_diff_this_difficulty_this_bin_cellBstimtype = ...
                        squeeze(mean(yyy_diff_hard_easy_cellBtB2Bk{j}(:,EI_t_range,dd,:),2));
                    EI_this_bin = yyy_diff_this_difficulty_this_bin_cellBstimtype(:,3)./...,
                        max(yyy_diff_this_difficulty_this_bin_cellBstimtype(:,1),...
                        yyy_diff_this_difficulty_this_bin_cellBstimtype(:,2));
                    
                    median_EI_this_bin = nanmedian(EI_this_bin);
                    mean_EI_this_bin = nanmean(EI_this_bin);
                    [~,p] = ttest(EI_this_bin,1);
                    
                    axes(h_subplot(ss));
                    hist(EI_this_bin,hist_range);    hold on; axis tight;
                    view(270,90); xlim(plot_range);
                    set(gca,'xaxislocation','top');
                    
                    if ss<sliding_n
                        %set(gca,'xticklabel',[],'yticklabel',[]);
                        set(gca,'xticklabel',[]);
                    end
                    
                    ylabel(sprintf('%g ~ %g',ttt(find(EI_t_range,1)),ttt(find(EI_t_range,1,'last'))));
                    title(sprintf('med %.3g\nmean %.3g\np = %.2g',median_EI_this_bin,mean_EI_this_bin,p));
                    
                    maxmaxy = max(maxmaxy,max(ylim));
                    
                    plot([median_EI_this_bin median_EI_this_bin],[0 20],'r-','linew',2);
                    plot([mean_EI_this_bin mean_EI_this_bin],[0 20],'g-','linew',2);
                    
                    % Save EI
                    EI_all_hardeasy_cellBtimeB2(:,ss,dd) = EI_this_bin;
                    
                end
                
                for ss = 1:sliding_n
                    %                     ylim(h_subplot(ss),[0 maxmaxy]);
                end
                
                SetFigure(15);
            end
            
            % -- Correlation between hard and easy EIs ---
            set(figure(750+dd+1),'unit','pix','name',['Enhance index, [hard, easy] = ' num2str(dd) ', ' methods_of_select{ms,2}],'pos',[2 49 1250 911]); clf
            [~,h_subplot] = tight_subplot(fix(sqrt(sliding_n)),ceil(sliding_n/fix(sqrt(sliding_n))),[0.05,0.02],[0.05,0.05]);
            
            for ss = 1:sliding_n % Sliding window
                
                % Easy vs Hard
                EI_t_range = (ss-1)*sliding_step + ttt_1 <= ttt & ttt < (ss-1)*sliding_step + ttt_1 + sliding_width;
                xx = EI_all_hardeasy_cellBtimeB2(:,ss,1);
                yy = EI_all_hardeasy_cellBtimeB2(:,ss,2);
                
                [~,p] = ttest(xx,yy);
                
                valid_EI = xx >0 & xx < 5 & yy >0 & yy < 5; % Just work-around, should select more carefully. @HH20160908
                h = LinearCorrelation(...
                    { xx(valid_EI), xx(~valid_EI)},...
                    { yy(valid_EI), yy(~valid_EI)},...
                    'Axes',h_subplot(ss),'Xlabel',sprintf('Hard (%g ~ %g), p= %g, x-y= %g',ttt(find(EI_t_range,1)),ttt(find(EI_t_range,1,'last')),p,mean(xx)-mean(yy)),...
                    'Ylabel','Easy','MarkerSize',5,'FaceColors',{'k','none'});
                delete([h.group(2).line]);
                axis([-1 5 -1 5]); axis square; legend off;
                hold on;
                plot(xlim,ylim,'--k');
                
                % Mean + sem
                mean_xx = mean(xx);
                sem_xx = std(xx)/sqrt(length(xx));
                mean_yy = mean(yy);
                sem_yy = std(yy)/sqrt(length(yy));
                plot([mean_xx-sem_xx mean_xx+sem_xx],[mean_yy mean_yy],'r-');
                plot([mean_xx mean_xx],[mean_yy-sem_yy mean_yy+sem_yy],'r-');
                
                % Show individual cell selected from the figure. HH20150424
                h_line = plot(xx,yy,'visible','off'); hold on;
                set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, methods_of_select{ms,1}});
                
                % EZZoom
                uicontrol('Style','pushbutton','String','tight','unit','norm','Position',[0.024 0.944 0.067 0.034],...
                    'callback','axis(findall(gcf,''type'',''axes''),''tight'')');
                uicontrol('Style','pushbutton','String','[-1 5]','unit','norm','Position',[0.096 0.944 0.067 0.034],...
                    'callback','axis(findall(gcf,''type'',''axes''),[-1 5 -1 5])');
                
            end
        end
    end

tuning_pack = [];

    function f1p2p7(debug)      % Rate 2.7. Tuning curve analysis. (Dora's tuning) 20161019
        if debug
            dbstack;
            keyboard;
        end
        
        %% Tuning: At 3 different phases. HH20150414
        % Now this part is dirty and quick.
        % Here are several things to be done:
        %  1. Align to stimulus onset AND saccade onset  (done)
        %  2. Deal with different heading angles
        %  3. Put this part elsewhere
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        j = 1;
        select_for_tuning = select_tcells;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % find_for_tuning = find(select_for_tuning);
        
        any_cell = any(group_ChoicePreference_pvalue(:,:,1) < 0.01,1)';
        
        %         any_cell = select_for_tuning;
        find_for_tuning = find(any_cell);
        
        
        for k = 1:3
                        tuning_cell = (group_ChoicePreference_pvalue(k,:,1) < 0.01)';
%             tuning_cell = (Choice_pref_p_value_all(k,:,1) < 0.01)';
            tmp = tuning_cell + select_for_tuning;
%             tmp(tmp == 0) = [];
            tuning_cell_in_any_cell{k} = find(tmp == 2);
            tuning_cell_index = ismember(find_for_tuning, tuning_cell_in_any_cell{k});
            tuning_cell_in_any_cell_index{k} = find(tuning_cell_index);
        end
        
        
        %         tuning_time_center = [605 750 1100 1500];
        tuning_time_center = [700 1000 1050 1500];
%         tuning_time_center = [735 925 1100 1500];
%                 tuning_time_center = [600 800 1000 1200];
        
        for tt = 1:length(tuning_time_center)
            [~,tmp_t] = min(abs(CP_ts{j} - tuning_time_center(tt)));
            tuning_time_phase(tt) = tmp_t;
        end
        %         tuning_time_phase_title = {'Vest largest', 'Comb largest','Vis largest', 'Pre-sac'};
        tuning_time_phase_title = {'Acceleration peak', 'speed peak','Deceleration peak', 'Stim end'};
        
        
        % Suppose all the unique_headings are the same!!
        unique_heading = group_result(representative_cell).mat_raw_PSTH.CP{j,1}.raw_CP_result{1}.Neu_tuning(:,1);
        % HH20160214: Patch calculation for LEFT & RIGHT choices of 0 headings
        zero_index = find(unique_heading == 0);
        unique_heading_two_zeros = [unique_heading(1:zero_index-1); -eps; eps ;unique_heading(zero_index+1:end)];
        %         unique_heading_two_zeros = [unique_heading(1:zero_index-1) ;unique_heading(zero_index+1:end)];  % No zero version
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parameters for tuning curve normalization
        
        %         time_for_dynamic_range = find(CP_ts{j} < time_markers{j}(1,2)); % all pre-sac times
        time_for_dynamic_range = [tuning_time_phase(1) - 5 tuning_time_phase(1) + 5]; % Around stimulus center
        modalities_share_same_dynamic_range = 1;
        linear_or_sigmoid = 1; % 1: linear fitting; 2: sigmoid fitting
        min_trials_dora_tuning = 3;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if isempty(tuning_pack) % Only run once if needed
            
            progressbar()
            % Pack tuning data
            for i = 1:sum(any_cell)  % For cells
                
                progressbar(i/sum(any_cell))
                
                % Only include three modalities
                if length(group_result(find_for_tuning(i)).mat_raw_PSTH.unique_stim_type) < 3
                    continue;
                end
                
                % Align preferred direcions to 'rightward'. So if PREF == 1 for a certain cell, we
                % flip all its tuning curve.
                flip_tuning = group_result(find_for_tuning(i)).PREF_PSTH == 1;
                
                % Find dynamic range
                dynamic_min_own_range = inf * [1 1 1]; % Each modality have its own dynamic range
                dynamic_max_own_range = -inf * [1 1 1];
                
                dynamic_min = inf; % Each modality have its own dynamic range
                dynamic_max = -inf;
                
                % Pack data (all tuning phases)
                
                
                for pp = 1:length(CP_ts{j})
                    
                    
                    for k = 1:3
                        this_raw = group_result(find_for_tuning(i)).mat_raw_PSTH.CP{j,k}.raw_CP_result{pp};
                        this_tuning = this_raw.Neu_tuning(:,2);
                        this_tuning_correctonly = this_raw.Neu_tuning_correctonly(:,2);
                        this_tuning_dora = this_raw.Neu_tuning_Dora_matrix_mean';  % HH20161019
                        this_tuning_dora_n = this_raw.Neu_tuning_Dora_matrix_n';
                        
                        % @HH20160214: Patch calculation for LEFT & RIGHT choices of 0 headings (because I failed to do this
                        % in the original CP_HH and I feel hesitant to redo all the batch files from A to Z right now...).
                        % Note, sadly, that this could be NaN because I also failed to pack the spike counts into
                        % spike_counts_allheadings_grouped if the number of left/right choices were fewer than 3... WTF...
                        
                        % @HH20160907. Today I redo all the batch files finally!
                        
                        % Discard dora tuning which has less than "min_trials_dora_tuning" trials. HH20161109
                        this_tuning_dora(this_tuning_dora_n<min_trials_dora_tuning) = nan;
                        
                        if length(this_tuning) == length(unique_heading)
                            zero_spikes_left_right = [mean(this_raw.spike_counts_allheadings_grouped{zero_index, 1}); mean(this_raw.spike_counts_allheadings_grouped{zero_index, 2})];
                            this_tuning_correctonly = [this_tuning_correctonly(1:zero_index-1); zero_spikes_left_right; this_tuning_correctonly(zero_index+1:end)];
                        else % No zero heading at all in the file
                            zero_spikes_left_right = [nan;nan];
                            this_tuning_correctonly = [this_tuning_correctonly(1:zero_index-1); zero_spikes_left_right; this_tuning_correctonly(zero_index:end)];
                            this_tuning = [this_tuning(1:zero_index-1); nan ;this_tuning(zero_index:end)];
                            this_tuning_dora = [this_tuning_dora(1:zero_index-1,:); nan nan; this_tuning_dora(zero_index:end,:)];
                        end
                        
                        %                     this_tuning_correctonly = [this_tuning_correctonly(1:zero_index-1); this_tuning_correctonly(zero_index+1:end)]; % No zero version
                        
                        % Align preferred direcions to 'rightward'. So if PREF == 1 for a certain cell, we
                        % flip all its tuning curve.
                        if flip_tuning
                            this_tuning = flipud(this_tuning);
                            this_tuning_correctonly = flipud(this_tuning_correctonly);
                            
                            % For the dora tuning matrix, if PREF == left, should flip both lr and ud.HH20161019
                            %               -    0    +                             null   0    pref
                            %   L choice    xx   xx   nan     ==>   null choice      xx     xx    nan
                            %   R choice    nan  xx  xx             pref choice      nan    xx    xx
                            this_tuning_dora = flipud(this_tuning_dora);
                            this_tuning_dora = fliplr(this_tuning_dora);
                        end
                        
                        % Update dynamic range
                        if sum(pp == time_for_dynamic_range) > 0
                            dynamic_min = min(dynamic_min, min(this_tuning_correctonly(:)));
                            dynamic_max = max(dynamic_max, max(this_tuning_correctonly(:)));
                            
                            dynamic_min_own_range(k) = min(dynamic_min_own_range(k), min(this_tuning_correctonly(:)));
                            dynamic_max_own_range(k) = max(dynamic_min_own_range(k), max(this_tuning_correctonly(:)));
                        end
                        
                        % Pack data
                        try
                            tuning_pack{1}{k,pp}(:,i) = this_tuning';
                            tuning_pack{2}{k,pp}(:,i) = this_tuning_correctonly';
                            tuning_pack{3}{k,pp}(:,:,i) = this_tuning_dora';  % HH20161019
                        catch
                            tuning_pack{1}{k,pp}(:,i) = nan;
                            tuning_pack{2}{k,pp}(:,i) = nan;
                            tuning_pack{3}{k,pp}(:,:,i) = nan;
                            keyboard;
                        end
                    end
                end
                
                % Normalization according to dynamic range of each neuron
                % (lowest in all tuning curve = 0, highest = 1)
                offset = dynamic_min;
                gain = dynamic_max - dynamic_min;
                
                offset_own_range = dynamic_min_own_range;
                gain_own_range = dynamic_max_own_range - dynamic_min_own_range;
                
                for pp = 1:length(CP_ts{j})
                    for k = 1:3
                        if modalities_share_same_dynamic_range
                            % Three modalities share the same dynamic range
                            tuning_pack{1}{k,pp}(:,i) = tuning_pack{1}{k,pp}(:,i) - offset;
                            tuning_pack{1}{k,pp}(:,i) = tuning_pack{1}{k,pp}(:,i) / gain;
                            
                            tuning_pack{2}{k,pp}(:,i) = tuning_pack{2}{k,pp}(:,i) - offset;
                            tuning_pack{2}{k,pp}(:,i) = tuning_pack{2}{k,pp}(:,i) / gain;
                            
                            tuning_pack{3}{k,pp}(:,:,i) = tuning_pack{3}{k,pp}(:,:,i) - offset;
                            tuning_pack{3}{k,pp}(:,:,i) = tuning_pack{3}{k,pp}(:,:,i) / gain;  % HH20161019
                        else
                            % Three modalities have their own dynamic ranges
                            tuning_pack{1}{k,pp}(:,i) = tuning_pack{1}{k,pp}(:,i) - offset_own_range(k);
                            tuning_pack{1}{k,pp}(:,i) = tuning_pack{1}{k,pp}(:,i) / gain_own_range(k);
                            
                            tuning_pack{2}{k,pp}(:,i) = tuning_pack{2}{k,pp}(:,i) - offset_own_range(k);
                            tuning_pack{2}{k,pp}(:,i) = tuning_pack{2}{k,pp}(:,i) / gain_own_range(k);
                            
                            tuning_pack{3}{k,pp}(:,:,i) = tuning_pack{3}{k,pp}(:,:,i) - offset_own_range(k);  % HH20161019
                            tuning_pack{3}{k,pp}(:,:,i) = tuning_pack{3}{k,pp}(:,:,i) / gain_own_range(k);
                        end
                    end
                end
                
            end
        end
        
        % Get means and sems
        tuning_mean_all = nan(length(unique_heading),length(CP_ts{j}),3);
        tuning_sem_all = tuning_mean_all;
        tuning_mean_correctonly = nan(length(unique_heading_two_zeros),length(CP_ts{j}),3);
        tuning_sem_correctonly = tuning_mean_correctonly;
        tuning_mean_dora = nan(2,length(unique_heading),length(CP_ts{j}),3);
        tuning_sem_dora = tuning_mean_dora;
        tuning_n_dora = tuning_mean_dora;
        tuning_sig_fit_correctonly = nan(3,length(CP_ts{j}),3);
        %         tuning_linear_fit_correctonly = nan(4,length(CP_ts{j}),3,2);
        tuning_linear_fit_correctonly = nan(4,length(CP_ts{j}),3,3);  % Changed by ZZ @ 20230619
        % for the last dimension, the first one is for LEFT headings,
        % the second for Right headings, the third for abs headings
        
        
        % For sigmoid fitting. HH20150528
        %         sigfunc = @(sig_para, x)(sig_para(1)+ sig_para(2)./ (1 + exp(-x/sig_para(3))));
        sigfunc = @(sig_para, x)(sig_para(1)+ sig_para(2) * normcdf(x,0,sig_para(3)));
        
        function err = cost_function(q,data_cum)
            x = data_cum(:,1);
            y = data_cum(:,2);
            
            z = sigfunc(q,x);
            err = norm(z-y);
        end
        
        % Plotting tuning curves at three time points
        set(figure(145506),'name',['Dora tuning, j = ' num2str(j)]); clf;
        set(gcf,'uni','norm','pos',[0.001       0.078       0.835       0.838]);
        
        for pp = 1:length(CP_ts{j})
            
            for k = 1:3
                % Mean and sem
                % Patch calculation for LEFT & RIGHT choices of 0 headings. HH20160214
                this_tuning_all = tuning_pack{1}{k,pp}(:,~isnan(tuning_pack{1}{k,pp}(1,:)));
                tuning_mean_all(:,pp,k) = nanmean(this_tuning_all,2);
                tuning_sem_all(:,pp,k) = nanstd(this_tuning_all,[],2)./sqrt(sum(~isnan(this_tuning_all),2));
                
                % this_tuning_correctonly =  tuning_pack{2}{k,pp}(:,~isnan(tuning_pack{2}{k,pp}(1,:)));
                this_tuning_correctonly =  tuning_pack{2}{k,pp}(:,  tuning_cell_in_any_cell_index{k}); % HH20180612  % Only select significant ChoicePref cells in each conditions
                
                tuning_mean_correctonly(:,pp,k) = nanmean(this_tuning_correctonly,2);
                tuning_sem_correctonly(:,pp,k) = nanstd(this_tuning_correctonly,[],2)./sqrt(sum(~isnan(this_tuning_correctonly),2));
                
                this_tuning_dora =  tuning_pack{3}{k,pp}(:,:,~isnan(tuning_pack{2}{k,pp}(1,:))); % Use the same NaN indicator
                tuning_mean_dora(:,:,pp,k) = nanmean(this_tuning_dora,3);
                tuning_n_dora(:,:,pp,k) = sum(~isnan(this_tuning_dora),3);
                tuning_sem_dora(:,:,pp,k) = nanstd(this_tuning_dora,[],3)./sqrt(tuning_n_dora(:,:,pp,k));
                
                
                if linear_or_sigmoid == 2   % Fitting sigmoid function. HH20150528
                    yy = nanmean(this_tuning_correctonly,2);
                    %                 tuning_sig_fit_correctonly(:,pp,k) = nlinfit(unique_heading, yy, sigfunc, [min(yy) range(yy) 1]);
                    
                    % Begin optimization
                    quick = fminsearch(@(q)cost_function(q,[unique_heading_two_zeros yy]),[min(yy) range(yy) 1]);
                    
                    % Output
                    tuning_sig_fit_correctonly(:,pp,k) = quick;
                    
                elseif linear_or_sigmoid == 1    % Linear fitting using group data. HH20161109
                    hh_index = {1:zero_index zero_index+1:size(this_tuning_correctonly,1)};
                    
                    for hh = LEFT:RIGHT % LEFT, RIGHT
                        
                        yy = this_tuning_correctonly(hh_index{hh},:)';
                        xx = unique_heading(hh_index{hh}-(hh==RIGHT))';
                        xx = repmat(xx,size(yy,1),1);
                        xx = xx(~isnan(yy)); yy = yy(~isnan(yy));
                        
                        [r,p] = corr(xx,yy,'type','pearson');
                        [para,~] = polyfit(xx,yy,1);
                        
                        tuning_linear_fit_correctonly(:,pp,k,hh) = [r p para];
                        %                         tuning_linear_fit_correctonly([3 4],pp,k,hh) = para;
                        
                    end
                    
                    % I use the "delta norm firing as a function of |heading|"
                    % to combine the two sides. It works!   HH20180606
                    this_tuning_diff = this_tuning_correctonly(hh_index{2},:) - ...
                        flipud(this_tuning_correctonly(hh_index{1},:));
                    yy_diff = this_tuning_diff';
                    xx_abs = repmat(unique_heading(hh_index{2}-1)',size(yy_diff,1),1);
                    xx_abs = xx_abs(~isnan(yy_diff)); yy_diff = yy_diff(~isnan(yy_diff));
                    [r,p] = corr(xx_abs(:),yy_diff(:),'type','pearson');
                    
                    [para, ~] = polyfit(xx_abs(:), yy_diff(:),1);
                    tuning_linear_fit_correctonly(:,pp,k,3) = [r p para];
                    
                    
                    %                     % Share the same p
                    %                     tuning_linear_fit_correctonly([1 2],pp,k,LEFT) = [r, p]';
                    %                     tuning_linear_fit_correctonly([1 2],pp,k,RIGHT) = [r, p]';
                    
                end
            end
        end
        
        % ---  Plotting ---
        for pp = 1:length(tuning_time_phase)
            
            for k = 1:3  % For each stim type
                
                % Plotting
                subplot(3,length(tuning_time_phase),pp); hold on; ylabel('All');
                
                % Traditional all heading
                %{
                plot(unique_heading,tuning_mean_all(:,tuning_time_phase(pp),k),'o','markersize',9,'color',colors(k,:),'LineWid',2);
                h = errorbar(unique_heading,tuning_mean_all(:,tuning_time_phase(pp),k),tuning_sem_all(:,tuning_time_phase(pp),k),'color',colors(k,:),'LineWid',2);
                errorbar_tick(h,10000);
                %}
                
                % Dora tuning
                
                hh_index = {1:zero_index zero_index:length(unique_heading)};
                cc_marker = {'<','>'};
                
                for cc = 1:2
                    for hh = LEFT:RIGHT
                        
                        h = scatter(unique_heading(hh_index{hh})+0.2*sign(1.5-cc),...
                            tuning_mean_dora(cc,hh_index{hh},tuning_time_phase(pp),k),...
                            tuning_n_dora(cc,hh_index{hh},pp,k),...
                            condition_colors(k,:),[cc_marker{cc}]);
                        
                        if cc == hh
                            set(h,'markerfacecolor',condition_colors(k,:));
                        else
                            set(h,'markerfacecolor','none');
                        end
                    end
                    
                    plot(unique_heading,tuning_mean_dora(cc,:,tuning_time_phase(pp),k),'color',condition_colors(k,:),'LineWid',2)
                end
                
                % plot(unique_heading,tuning_mean_dora(:,:,tuning_time_phase(pp),k),'o','markersize',9,'color',colors(k,:),'LineWid',2);
                % h = errorbar([unique_heading'+0.2; unique_heading'-0.2],tuning_mean_dora(:,:,tuning_time_phase(pp),k),tuning_sem_dora(:,:,tuning_time_phase(pp),k),'color',colors(k,:),'LineWid',2,'linestyle','none');
                
                title([tuning_time_phase_title{pp} ]);
                axis tight; xlim(xlim*1.1);
                
                % --------- Tuning all ---------
                subplot(3,length(tuning_time_phase),pp + length(tuning_time_phase));  hold on; ylabel('Correct + wrong');
                
                errorbar(unique_heading,tuning_mean_all(:,tuning_time_phase(pp),k),...
                    tuning_sem_all(:,tuning_time_phase(pp),k),'o-','color',condition_colors(k,:),'LineWid',2,'capsize',0);
                
                
                % --------- Correct only with fitting ----------
                subplot(3,length(tuning_time_phase),pp + 2*length(tuning_time_phase));  hold on; ylabel('Correct only');
                
                plot(unique_heading_two_zeros,tuning_mean_correctonly(:,tuning_time_phase(pp),k),'o','markersize',9,'color',condition_colors(k,:),'LineWid',2);
                
                errorbar(unique_heading_two_zeros,tuning_mean_correctonly(:,tuning_time_phase(pp),k),...
                    tuning_sem_correctonly(:,tuning_time_phase(pp),k),'color',condition_colors(k,:),'LineWid',2,'linestyle','none','capsize',0);
                
                % try errorbar_tick(h,10000); catch end;
                
                title([tuning_time_phase_title{pp} ', t = ' num2str(CP_ts{j}(tuning_time_phase(pp)))]);
                text(-8,0.4+0.05*k,sprintf('%g, %.2g', size(tuning_cell_in_any_cell{k},1), ...
                    tuning_linear_fit_correctonly(2,tuning_time_phase(pp),k,1)),'color',condition_colors(k,:));
                text(2,0.1+0.05*k,sprintf('%g, %.2g', size(tuning_cell_in_any_cell{k},1), ...
                    tuning_linear_fit_correctonly(2,tuning_time_phase(pp),k,2)),'color',condition_colors(k,:));
                
                if linear_or_sigmoid == 2        % Plot sigmoid fitting
                    xx = linspace(min(unique_heading),max(unique_heading),100);
                    plot(xx,sigfunc(tuning_sig_fit_correctonly(:,tuning_time_phase(pp),k),xx),'color',condition_colors(k,:),'linew',3);
                    
                elseif linear_or_sigmoid == 1   % Plot linear fitting
                    
                    line_type = {'-','--'}; % Significant/Insignificant
                    for hh = LEFT:RIGHT % LEFT, RIGHT
                        xx = unique_heading(hh_index{hh})';
                        para = tuning_linear_fit_correctonly(:,tuning_time_phase(pp),k,hh); % [r,p,b,k]
                        plot(xx,para(4)+para(3)*xx,'-','color',condition_colors(k,:),'linew',3,'linestyle',line_type{2-(para(2)<0.05)});
                    end
                end
            end
            
            axis tight; xlim(xlim*1.1);
            SetFigure(15);
            
        end
        
        % Linear fitting parameters
        set(figure(145525),'name',['Linear fittings, j = ' num2str(j)]); clf;  hold on
        set(gcf,'uni','norm','pos',[0.002       0.532       0.411        0.38]);
        
        for hh = 1:3
            subplot(3,1,hh); hold on;
            for k = 1:3
                
                %             hh = 2;
                plot(CP_ts{j}, tuning_linear_fit_correctonly(3,:,k,hh),'linew',3,'linestyle',line_type{1},'color',condition_colors(k,:));
                sig = tuning_linear_fit_correctonly(2,:,k,hh) < 0.01;
                plot(CP_ts{j}(sig), mean(tuning_linear_fit_correctonly(3,sig,k,hh),3),'o','linew',2,...
                    'markersize',7,'color',condition_colors(k,:));
                %                 plot(CP_ts{j}, mean(tuning_linear_fit_correctonly(2,:,k,hh),3),[line_type{3-hh}],...
                %                     'linew',1,'color',colors(k,:));
            end
            ylim([-0.03 0.05])
        end
        ylabel('Slope of linear fitting');
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{j}(1), 0  + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
        SetFigure(15);
        
    end


dora_tuning_mean_each_cell = []; % Cache the dora tuning for each cell calculated in f1p2p8
dora_tuning_sem_each_cell = []; % Cache the dora tuning for each cell calculated in f1p2p8
dora_tuning_n_each_cell = [];
partial_corr_timewins = {[0 1500],'0 ~ 1500 ms (all stim)';
    %                              [-200 0],'-300 ~ 0 ms (before stim)';
    [0 300],'0 ~ 300 ms (stim start)';
    [300 600],'300 ~ 600 (motion initiation)';
    [600 800],'600 ~ 800 ms (acc peak)'
    [800 1000],'800 ~ 1000 ms (vel peak)';
    [1000 1200], '1000 ~ 1200 ms (deAcc peak)';
    [1200 1500],'1200 ~ 1500 ms (stim end)'
    %                              [1500 1800],'1500 ~ 2000 ms (delay)';
    %                              [2000 2300],'2000 ~ 2300 ms (postsac)'
    };
select_for_partial = [];

    function f1p2p8(debug)    % Rate 2.8. Partial correlation analysis. (Dora's method) 20170327 @UNIGE
        if debug
            dbstack;
            keyboard;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        j = 1;
%         select_for_partial = select_tcells;
        select_for_partial = select_all;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        find_for_partial = find(select_for_partial);
        unique_heading = group_result(representative_cell).unique_heading;
        partial_corr_coef_all = nan(sum(select_for_partial),size(partial_corr_timewins,1),2,3);  % [cell number, time epoch, (heading coeff, choice coeff), stim_type]
        partial_corr_p_all = nan(sum(select_for_partial),size(partial_corr_timewins,1),2,3);
        partial_corr_coef_all_flipped = nan(sum(select_for_partial),size(partial_corr_timewins,1),2,3);  % For plotting partial corr over time
        
        anova2_p_all = nan(sum(select_for_partial),size(partial_corr_timewins,1),2,3);
        dora_tuning_mean_each_cell = nan(sum(select_for_partial),size(partial_corr_timewins,1),3,2,length(unique_heading));
        dora_tuning_sem_each_cell = nan(sum(select_for_partial),size(partial_corr_timewins,1),3,2,length(unique_heading));
        dora_tuning_n_each_cell = zeros(sum(select_for_partial),size(partial_corr_timewins,1),3,2,length(unique_heading));
        
        % Calculate partial correlation coefficients and p-values for each cell
        progressbar('cell num');
        for i = 1:sum(select_for_partial)  % For cells
            
            this_raw_spike_in_bin = group_result(find_for_partial(i)).mat_raw_PSTH.spike_aligned{1,j};
            this_time = group_result(find_for_partial(i)).mat_raw_PSTH.spike_aligned{2,j};
            this_stim_type = group_result(find_for_partial(i)).mat_raw_PSTH.stim_type_per_trial;
            this_heading = group_result(find_for_partial(i)).mat_raw_PSTH.heading_per_trial;
            this_choice = group_result(find_for_partial(i)).mat_raw_PSTH.choice_per_trial;
            
            for tt = 1:size(partial_corr_timewins,1)
                count_win = partial_corr_timewins{tt,1}(1) <= this_time & this_time <= partial_corr_timewins{tt,1}(2);
                for k = 1:3
                    if isempty(find(this_stim_type==k, 1)); continue; end
                    
                    % --- Partial correlation ---
                    X=[];
                    X(:,1) = sum(this_raw_spike_in_bin(this_stim_type==k,count_win),2)...
                        /range(partial_corr_timewins{tt,1})*1e3; % Average firing rate in Hz
                    X(:,2) = this_heading(this_stim_type==k);
                    X(:,3) = this_choice(this_stim_type==k);
                    
                    [r,p] = partialcorr(X);
                    partial_corr_coef_all(i,tt,:,k) = r(1,2:3);
                    partial_corr_p_all(i,tt,:,k) = p(1,2:3);
                    
                    % --- Anova2 p values (does not exclude interactions of heading and choice, but Dora used this) ---
                    anova2_p_all(i,tt,:,k) = anovan(X(:,1),{X(:,2),X(:,3)},'display','off')';
                    
                    % --- Dora tuning (putting it here is more flexbile than that in CP_HH calculation) ---
                    real_unique_this_heading = unique(this_heading);
                    withzero_unique_this_heading = unique([this_heading 0]);
                    for hh = 1:length(real_unique_this_heading)
                        for cc = LEFT:RIGHT
                            this_select = (X(:,2) == real_unique_this_heading(hh))&(X(:,3)==cc);
                            this_mean = mean(X(this_select,1));
                            this_sem = std(X(this_select,1))/sqrt(sum(this_select));
                            hh_shouldbe = withzero_unique_this_heading==real_unique_this_heading(hh); % Deal with cases without zero heading
                            dora_tuning_mean_each_cell(i,tt,k,cc,hh_shouldbe) = this_mean;
                            dora_tuning_sem_each_cell(i,tt,k,cc,hh_shouldbe) = this_sem;
                            dora_tuning_n_each_cell(i,tt,k,cc,hh_shouldbe) = sum(this_select);
                        end
                    end
                end
                
            end
            progressbar(i/sum(select_for_partial));
        end
        
        % Use R^2 (Fig.4 of Zaidel 2017)
        partial_corr_coef_all_flipped = partial_corr_coef_all.^2;
        
        
        %% Drawing
        set(figure(3099+figN),'name',sprintf('Partial correlation, j = %g, "any significant" out of N = %g, "%s" cells',...
            j,sum(select_for_partial),t_criterion_txt)); clf; figN = figN+1;
        set(gcf,'uni','norm','pos',[0       0.038       0.994       0.877]);
        [h_sub,~] = tight_subplot(3,size(partial_corr_timewins,1),[0.05 0.02]);
        
        for tt = 1:size(partial_corr_timewins,1)
            for k = 1:3
                
                axes(h_sub(k+(tt-1)*3));
                % Use partial corr p value
                %                 heading_sig = partial_corr_p_all(:,tt,1,k)<0.05;
                %                 choice_sig = partial_corr_p_all(:,tt,2,k)<0.05;
                
                % Use ANOVA2 p value (Dora)
                heading_sig = anova2_p_all(:,tt,1,k)<0.05;
                choice_sig = anova2_p_all(:,tt,2,k)<0.05;
                
                h1=plot(partial_corr_coef_all((heading_sig|choice_sig) & ~(heading_sig & choice_sig),tt,1,k),...
                    partial_corr_coef_all((heading_sig|choice_sig) & ~(heading_sig & choice_sig),tt,2,k),...  % "One sig" cells
                    'o','color',condition_colors(k,:));
                hold on;
                h2=plot(partial_corr_coef_all(heading_sig&choice_sig,tt,1,k),...
                    partial_corr_coef_all(heading_sig&choice_sig,tt,2,k),...  % "All sig" cells
                    'o','color','k','markerfacecol','k');
                
                xx = partial_corr_coef_all(heading_sig|choice_sig,tt,1,k);
                yy = partial_corr_coef_all(heading_sig|choice_sig,tt,2,k);
                h_all = plot(xx,yy,'visible','off');
                
                % Draw line if significant
                [rrr,ppp] = corr(xx,yy,'type','Pearson');
                if ppp < inf % 0.05
                    coeff = pca([xx yy]);
                    linPara(1) = coeff(2) / coeff(1);
                    linPara(2) = mean(yy)- linPara(1) *mean(xx);
                    
                    % -- Plotting
                    xxx = linspace(min(xx),max(xx),150);
                    yyy = linPara(1) * xxx + linPara(2);
                    plot(xxx,yyy,'linew',2,'color',condition_colors(k,:));
                end
                
                select_actual_plot = zeros(N,1);
                find_for_partial = find(select_for_partial);
                select_actual_plot(find_for_partial(heading_sig|choice_sig)) = 1;
                set([gca h_all],'ButtonDownFcn',{@Show_individual_cell, h_all, select_actual_plot});
                
                
                axis([-1 1 -1 1]);
                hold on; plot([-1 1],[0 0],'k--'); plot([0 0],[-1 1],'k--');
                text(-0.9,-0.8,sprintf('%g,%g,%g\nr = %.3g\np= %.3g',sum(heading_sig),sum(choice_sig),sum(heading_sig&choice_sig),rrr,ppp));
                
                if k == 1
                    title(partial_corr_timewins{tt,2});
                end
                if tt == 1 && k == 1
                    xlabel('Partial corr (Heading)');
                    ylabel('Partial corr (Choice)');
                end
                
                
            end
        end
        
        %% Plot partial correlation over time (Zaidel 2017, Figure 6)  HH20180821
        set(figure(194959),'name',sprintf('Partial correlation over time, j = %g, "any significant" out of N = %g, "%s" cells',...
            j,sum(select_for_partial),t_criterion_txt)); clf; figN = figN+1;
        set(gcf,'uni','norm','pos',[0       0.038       0.994       0.877]);
        
        time_centers = mean(reshape([partial_corr_timewins{2:end,1}],2,[]),1);
        
        for k = 1:3
            subplot(1,3,k);
            for pp = 1:2
                
                aver_this = nanmean(partial_corr_coef_all_flipped(:,2:end,pp,k),1); % Note the first time epoch is the whole trial
                sem_this = nanstd(partial_corr_coef_all_flipped(:,2:end,pp,k),[],1)/sqrt(size(partial_corr_coef_all_flipped,1));
                
                if pp == 1
                    col = condition_colors(k,:);
                else
                    col = [0 0 0];
                end
                errorbar(time_centers, aver_this, sem_this,'color',col,'linew',2); hold on;
                
            end
        end
        title('Temporal resolution is too low here for Fisher info. I put the real one in f6p5r2.');
        
    end

    function f1p3(debug)      % Rate 3. Correct / Wrong Trials
        if debug
            dbstack;
            keyboard;
        end
        %%
        %         methods_of_select = {
        %             %select_bottom_line, 'All cells';
        %             select_tcells, 'Typical cells';
        %             % select_no_tcells, 'Non-typical cells'
        %             };
        methods_of_select = {
            % select_bottom_line, 'All cells';
            selected_t_vest, 'Vest. cell';
            selected_t_vis, 'Vest. cell';
            selected_t_comb, 'Vest. cell';
            select_tcells, 'Typical cells';
            };
        
        
        unique_heading = group_result(representative_cell).mat_raw_PSTH.CP{j,1}.raw_CP_result{1}.Neu_tuning(:,1);
        unique_heading_for_correct_wrong = unique_heading(unique_heading > 0);
        
        for ms = 1 %:size(methods_of_select,1)
            set(figure(2254-ms),'name',['Average PSTH (Correct vs Wrong), ' methods_of_select{ms,2}],'pos',[18 67 1645 898]); clf
            h_subplot = tight_subplot(3,1 + length(unique_heading_for_correct_wrong),[0.05 0.02],[0.1 0.1]);
            
            j = 1;
            for k = 1:3
                
                % ------ All correct and wrong ------ %
                h = SeriesComparison({PSTH_outcomes_Norm{1}(methods_of_select{k,1},:,:,k) PSTH_outcomes_Norm{2}(methods_of_select{k,1},:,:,k) PSTH_outcomes_Norm{3}(methods_of_select{k,1},:,:,k)},...
                    {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
                    'Colors',{condition_colors(k,:), condition_colors(k,:), 'k', 'k'},'LineStyles',{'-','--'},...
                    'ErrorBar',0,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(k));
                
                if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                
                % legend off;
                
%                 if k == 1
                    title([methods_of_select{k,2} ', n = ' num2str(sum(methods_of_select{k,1}))]);
%                 end
                
                %                 for tt = 1:3
                %                     plot([1 1] * time_markers{j}(1,tt),[0 1],'k','linestyle',marker_for_time_markers{j}{tt},'linew',2);
                %                 end
                
               xlim([rate_ts{j}(10) rate_ts{j}(end-10)]);  %ylim([.0 .7]);
                
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                
                legend off;
                
                % ------ Correct and wrong for each heading ------ %    @HH20150523
                
                for hh = 1:length(unique_heading_for_correct_wrong)
                    
                    % Reorganize to fit "correct +, correct -, wrong +, wrong -" for each heading
                    correct_pref_ind = 2 + hh*2-1 ; % Index in PSTH_correct_angles_norm (exclude first two 0 heading)
                    correct_null_ind = 2 + hh*2 ;  % Index in PSTH_correct_angles_norm
                    % Heading order in wrong trials changed in the main
                    % analysis code, so change the following part accordingly
                    wrong_pref_ind =  hh*2-1 ; % Index in PSTH_wrong_angles_norm
                    wrong_null_ind =  hh*2 ; % Index in PSTH_wrong_angles_norm
                    
                    for jjj = 1:3
                        PSTH_correct_wrong_this{jjj} = ...
                            cat(3, PSTH_correct_angles_Norm{jjj}(methods_of_select{k,1},:,[correct_pref_ind correct_null_ind],k),...
                            PSTH_wrong_angles_Norm{jjj}(methods_of_select{k,1},:,[wrong_pref_ind wrong_null_ind],k));
                    end
                    
                    SeriesComparison({PSTH_correct_wrong_this{1} PSTH_correct_wrong_this{2} PSTH_correct_wrong_this{3}},...
                        {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
                        'Colors',{condition_colors(k,:), condition_colors(k,:), 'k', 'k'},'LineStyles',{'-','--'},...
                        'ErrorBar',0,'Xlabel',[],'Ylabel',[],'axes',h_subplot(hh * 3 + k));
                   
                    
                    if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
                    %                     set(gca,'ytick',[]);
                    
                    legend off;
                    
%                     if k == 1
                        title([methods_of_select{k,2} ', ' num2str(unique_heading_for_correct_wrong(hh))  ' n = ' num2str(sum(methods_of_select{k,1}))]);
%                     end
                    
                   xlim([rate_ts{j}(10) rate_ts{j}(end-10)]);  % ylim([.0 .7]);
                    
                    % Gaussian vel
                    plot(Gauss_vel(:,1) + time_markers{j}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
                end
                
            end
            
            SetFigure(10);
        end
    end

marker_size = 15;
tcell_cross_size = 15;
tWin4CD = {   % Artificially choose 300-ms period for averaging ChoiceDiv for each condition according to divergence time
    [600 900], 'Vestibular period';
    [900 1200], 'Visual period';
    [750 1050], 'Combined period'
    };

    function f2p1(debug)      % ROC 1. CP, CDiv, MDiv
        if debug
            dbstack;
            keyboard;
        end
        
        %%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         select_for_CP = select_bottom_line;
        %         select_for_div = select_bottom_line;
        
        select_for_CP = select_tcells;
        select_for_ChoDiv = select_tcells;
        select_for_ModDiv = select_bottom_line;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% CP: Time-course
        
        set(figure(2099+figN),'name','CP','pos',[287 483 1161 475]); figN = figN+1; clf;
        
        for j = 1:3
            subplot(1,3,j);
            
            for k = 1:3
                
                ys = nanmean(CP{j}(select_for_CP,:,k));
                errors = nanstd(CP{j}(select_for_CP,:,k))/sqrt(sum(select_for_CP));  % change from mean and std to nanmean and nanstd
                h = shadedErrorBar(CP_ts{j},ys,errors,'lineprops',{'Color',condition_colors(k,:)},'transparent',transparent);
                set(h.mainLine,'LineWidth',2);
                hold on;
                
            end
            
            xlabel(['Center of ' num2str(group_result(representative_cell).mat_raw_PSTH.binSize_CP) ' ms time window']);
            
            title(sprintf('Grand CP (N = %g)',sum(select_for_CP)));
            xlim([CP_ts{j}(1) CP_ts{j}(end)]);
            ylim([0.4 0.7]);
            
            %             for tt = 1:3
            %                 plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',1.5);
            %             end
            if j == 1
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/5 + min(ylim),'--','linew',3,'color',[0.6 0.6 0.6]);
            end
            
            plot(xlim,[0.5 0.5],'k--');
            
        end
        
        SetFigure(17);
        
        %% ModDiv: Time-course
        ModDiv_linesty = {'-', '--', ':'};
        lgd = {'2-1'; '3-1'; '3-2'};
        
        SeriesComparison({ModDiv_All{1}(select_for_ModDiv,:,:) ModDiv_All{2}(select_for_ModDiv,:,:) ModDiv_All{3}(select_for_ModDiv,:,:)},...
            {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
            'Colors',{[0.5 0.5 0.5],[0.5 0.5 0.5],[0.5 0.5 0.5]},'LineStyles',ModDiv_linesty,...
            'ErrorBar',6,'Xlabel',[],'Ylabel',[],'transparent',transparent,...
            'CompareIndex',[1:3;1:3],...
            'CompareColor',[mat2cell(condition_colors,ones(3,1))],...
            'figN', figN,  'PCritical', 0.01);
        
        %         xlim([-300 2300]); ylim([-0.1 0.7 ]);
        legend off;  plot(xlim,[0 0],'k--');
        title(sprintf('Modality divergence (SU, N = %g)',sum(select_for_ModDiv)));
        set(gcf,'unit','norm','pos',[0.00952380952380952 0.521904761904762 0.315476190476191 0.385714285714286]);
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{1}(1),Gauss_vel(:,2)*range(ylim)/5 + min(ylim),'--','linew',2,'color',[0.6 0.6 0.6]);
        
        SetFigure(17);
        
        
        %% ChoiceDiv. HH20160213
        
        SeriesComparison({ChoiceDiv_All{1}(select_for_ChoDiv,:,:) ChoiceDiv_All{2}(select_for_ChoDiv,:,:) ChoiceDiv_All{3}(select_for_ChoDiv,:,:)},...
            {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
            'Colors',{condition_colors(1,:),condition_colors(2,:),condition_colors(3,:)},'LineStyles',{'-'},...
            'ErrorBar',6,'Xlabel',[],'Ylabel',[],'transparent',transparent,...
            'CompareIndex',[1:3;1:3],...
            'CompareColor',[mat2cell(condition_colors,ones(3,1))],...
            'PCritical', 0.01);
        
        %         xlim([-300 2300]); ylim([-0.1 0.7 ]);
        legend off;  plot(xlim,[0 0],'k--');
        title(sprintf('Choice divergence (SU, N = %g)',sum(select_for_ChoDiv)));
        set(gcf,'unit','norm','pos',[0.00952380952380952 0.521904761904762 0.315476190476191 0.385714285714286]);
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{1}(1),Gauss_vel(:,2)*range(ylim)/5 + min(ylim),'--','linew',2,'color',[0.6 0.6 0.6]);
        
        SetFigure(17);
        
        %% Linear fitting combined trace
        for j = 1:3
            % Choice divergence
            for k = 1:3
                selectCells_notNaN =  select_for_ChoDiv & (~isnan(ChoiceDiv_All{j}(:,1,k)));
                
                ys_CD{j}(k,:) = mean(ChoiceDiv_All{j}(selectCells_notNaN,:,k));
                errors = std(ChoiceDiv_All{j}(selectCells_notNaN,:,k))/sqrt(sum(selectCells_notNaN));
            end
        end
        
        ts = rate_ts{1};
        t_select = rate_ts{1}> 0 & rate_ts{1} <=1500;
        
        ramping1 = ys_CD{1}(1,t_select);
        ramping2 = ys_CD{1}(2,t_select);
        ramping3 = ys_CD{1}(3,t_select);
        
        w = fminsearch(@(w) sum((w(1)*ramping1 + w(2)*ramping2 - ramping3).^2), [.5 .5]);
        
        figure(2308); clf
        plot(ts(t_select),ramping1,'b',ts(t_select),ramping2,'r',...
            ts(t_select),ramping3,'g',ts(t_select),ramping1*w(1)+ramping2*w(2),'k');
        title(sprintf('Fitting Choice Div: w1 = %g, w2 = %g',w(1),w(2)))
        
        
    end


    function f2p2(debug)      % ROC 2. CDiv: Multisensory Enhancement
        if debug
            dbstack;
            keyboard;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         select_for_div =  select_bottom_line;
        select_for_div =  select_tcells;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% CDiv: Multisensory Enhancement
        % select_bottom_line = ([group_result.repN]' >= 8) & (xls_num{1}(:,header.Chan1) > 0) & (xls_num{1}(:,header.HD_TargFirst)~=0);% ...
        % & (~isnan(ChoiceDiv_All(:,1,k))) ;% & MemSac_DDI(:,4)<0.55 ;
        
        set(figure(2099+figN),'name','Multisensory Enhancement of CDiv'); clf; figN = figN+1;
        set(gcf,'uni','norm','pos',[0.044        0.49       0.604       0.289]);
        
%         ind = rate_ts{1} >= -175 & rate_ts{1} <= 1675; 
        
        for k = 1:3
            h = subplot(1,3,k);
            
            select_this_k = select_for_div & (~isnan(ChoiceDiv_ModDiffer{1}(:,1,k)));
            
            SeriesComparison({ChoiceDiv_ModDiffer{1}(select_this_k,:,k) ChoiceDiv_ModDiffer{2}(select_this_k,:,k) ChoiceDiv_ModDiffer{3}(select_this_k,:,k)},...
                {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
                'ErrorBar',6,...
                'CompareIndex',[1;1],'CompareColor',{condition_colors(k,:)},'PlotPs',0,'PCritical',0.01,...
                'Colors',{condition_colors(k,:)},'YLim',[-0.3 0.45],...
                'Transparent',transparent,'LineStyles',{'-'},'axes',h);
            
            legend off; xlim([-100 1600]);
            plot(xlim,[0 0],'k--');
            
            % Gaussian vel
            plot(Gauss_vel(:,1) ,Gauss_vel(:,2)*range(ylim)/5 ,'--','linew',2,'color',[0.6 0.6 0.6]);
            
        end
        
        title(sprintf('3-1 3-2 1-2 (N = %g)',sum(select_this_k)));
        SetFigure(15);
    end

% Global and updating variable the following part will use 
        ChoiceDiv_Error;
        ChoiceDiv_Correct;
        ChoiceDiv_CorrectMinusError;
        ChoiceDiv_Difficult;
        ChoiceDiv_Easy;
        ChoiceDiv_EasyMinusDifficult; 
        ChoicePref_Error; 
        ChoicePref_Difficult;
        ChoicePref_Easy;
        ChoicePref_Correct;

    function f2p3(debug)      % ROC 3. Easy and Difficult
        if debug
            dbstack;
            keyboard;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         select_for_div =  select_tcells & select_bottom_line;
for k = 1:3 % conditions
%         cpref_sig = any(Choice_pref_p_value_all(:,:,1) < 0.01); 
        cpref_sig{k} = Choice_pref_p_value_all(k,:,1) < 0.01; 

end
        select_for_div = cpref_sig; 
        %         select_for_div =  select_bottom_line;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% CDiv: Easy AND difficult
        set(figure(figN),'pos',[12 276 1272 684],'name','Easy & difficult of CDiv'); clf; hold on;
        
        % Concatenate Correct and error trials for plotting
        ChoiceDiv_EasyVSDiffic = cellfun(@(x,y) cat(3, x, y), ChoiceDiv_Easy, ChoiceDiv_Difficult, 'UniformOutput',0);
        SeriesComparison({ChoiceDiv_EasyVSDiffic{1}(select_for_div{1}',:,:) ChoiceDiv_EasyVSDiffic{2}(select_for_div{2}',:,:) ChoiceDiv_EasyVSDiffic{3}(select_for_div{3}',:,:)},...
            {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
            'Colors',mat2cell(cat(1, condition_colors, condition_colors), ones(6,1)),...
            'LineStyles', {'-','-','-','--','--','--'},...
            'ErrorBar',6, 'Xlabel',[], 'Ylabel', 'Choice Divergence', 'figN', figN,...
            'CompareIndex',[1,2,3; 4,5,6], 'CompareColor', mat2cell(condition_colors, ones(3,1)),...
            'Transparent', transparent, 'PCritical', 0.01); figN=figN+1;
        
        legend off;
        plot(xlim,[0 0],'k--');
        xlim([-200 1500+50])
%         title(['n= ' num2str(sum(select_for_div))]);
        
        SetFigure(20);
        
        %% CDiv: Easy-difficult
        
        % select_bottom_line = ([group_result.repN]' >= 8) & (xls_num{1}(:,header.Chan1) > 0) & (xls_num{1}(:,header.HD_TargFirst)~=0);% ...
        % & (~isnan(ChoiceDiv_All(:,1,k))) ;% & MemSac_DDI(:,4)<0.55;
        
        set(figure(figN),'name','Easy-difficult of CDiv'); clf; hold on;
        
        SeriesComparison({ChoiceDiv_EasyMinusDifficult{1}(select_for_div{1}',:,:) ChoiceDiv_EasyMinusDifficult{2}(select_for_div{2}',:,:) ChoiceDiv_EasyMinusDifficult{3}(select_for_div{3}',:,:)},...
            {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
            'Xlabel',[], 'Ylabel', 'Delta Choice Divergence',...
            'Colors',mat2cell(condition_colors,ones(3,1)),'ErrorBar',6,'PCritical',0.01,...
            'CompareIndex',[1 2 3; 1 2 3],'CompareColor',mat2cell(condition_colors,ones(3,1)),...
            'Transparent',transparent,'LineStyles',{'-','-','-'}, 'figN', figN);  figN = figN+1;
        
        plot(Gauss_vel(:,1), Gauss_vel(:,2)*range(ylim)/5 + 0,'--','linew',2,'color',[0.6 0.6 0.6]);
        
        legend off;
        plot(xlim,[0 0],'k--');
        xlim([-200 1500+50])

%         title(['n= ' num2str(sum(select_for_div))]);
        
        SetFigure(20);
        
    end

    function f2p4(debug)      % ROC 3. Correct and error
        if debug
            dbstack;
            keyboard;
        end
        
        %         select_for_div =  select_tcells;
        for k =1:3  % condition
            cpref_sig{k} = Choice_pref_p_value_all(k,:,1) < 0.01;

%             cpref_sig = any(Choice_pref_p_value_all(:,:,1) < 0.01);
        end
        select_for_div = cpref_sig;
        % select_for_div = {selected_t_vest; selected_t_vis;selected_t_comb};
        
        set(figure(figN),'pos',[12 276 1272 684],'name','Correct & Error of CDiv'); clf; hold on;
        
        % Concatenate Correct and error trials for plotting
        ChoiceDiv_CorrectVSError = cellfun(@(x,y) cat(3, x, y), ChoiceDiv_Correct, ChoiceDiv_Error, 'UniformOutput',0);
        SeriesComparison({ChoiceDiv_CorrectVSError{1}(select_for_div{1}',:,:) ChoiceDiv_CorrectVSError{2}(select_for_div{2}',:,:) ChoiceDiv_CorrectVSError{3}(select_for_div{3}',:,:)},...
            {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
            'Colors',mat2cell(cat(1, condition_colors, condition_colors), ones(6,1)),...
            'LineStyles', {'-','-','-','--','--','--'},...
            'ErrorBar',6, 'Xlabel',[], 'Ylabel', 'Choice Divergence', 'figN', figN,...
            'CompareIndex',[1,2,3; 4,5,6], 'CompareColor', mat2cell(condition_colors, ones(3,1)),...
            'Transparent', transparent, 'PCritical', 0.01); figN=figN+1;
        
        legend off;
        plot(xlim,[0 0],'k--');
%         title(['n= ' num2str(sum(select_for_div))]);
        xlim([-200 1500+50])
        SetFigure(20);
        
        %% ChoiceDiv Correct - Error
        set(figure(figN),'name','Correct-Error of CDiv'); clf; hold on;
        
        SeriesComparison({ChoiceDiv_CorrectMinusError{1}(select_for_div{1}',:,:) ChoiceDiv_CorrectMinusError{2}(select_for_div{2}',:,:) ChoiceDiv_CorrectMinusError{3}(select_for_div{3}',:,:)},...
            {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
            'Xlabel',[], 'Ylabel', 'Delta Choice Divergence',...
            'Colors',mat2cell(condition_colors,ones(3,1)),'ErrorBar',6,'PCritical',0.01,...
            'CompareIndex',[1 2 3; 1 2 3],'CompareColor',mat2cell(condition_colors,ones(3,1)),...
            'Transparent',transparent,'LineStyles',{'-','-','-'}, 'figN', figN);figN=figN+1;
        
        plot(Gauss_vel(:,1), Gauss_vel(:,2)*range(ylim)/5 + 0,'--','linew',2,'color',[0.6 0.6 0.6]);
        legend off;
        plot(xlim,[0 0],'k--');
        xlim([-200 1500+50])
%         title(['n= ' num2str(sum(select_for_div))]);
        
        SetFigure(20);
        
    end

    function f3p2p4(debug)    % Choice Preference: Correct and Error
        if debug
            dbstack;
            keyboard;
        end
        
                j = 1; % Stimulus period
        
        select_for_choicepref =  (squeeze(Choice_pref_p_value_all(:,:,j)) < 0.01)';
        
        monkeys = xls_num{1}(:,header.Monkey);  % HH20160613
        monkey1 = monkeys == 15; monkey1 = monkey1(select_bottom_line);
        monkey2 = monkeys == 13; monkey2 = monkey2(select_bottom_line);
        
%         set(figure(figN),'name','ChoicePref Correct vs. Error'); clf; hold on; axis off; 
%         ha = axes('position', [0.1 0.1 0.8 0.8]); 
        
        for st = 1:3
            
            set(figure(figN),'name','ChoicePref Correct vs. Error'); clf; hold on; axis off;
            ha = axes('position', [0.1 0.1 0.8 0.8]);
            
            % One-side paired t-test
            % Because I assume the ChoicePref in correct trials is larger
            % so I use one-side ttest to add more information 
            [~,p1] = ttest(abs(ChoicePref_Error(select_for_choicepref(:,st)&monkey1,j,st)), abs(ChoicePref_Correct(select_for_choicepref(:,st)&monkey1,j,st)) , 'Tail', 'left');
            [~,p2] = ttest(abs(ChoicePref_Error(select_for_choicepref(:,st)&monkey2,j,st)), abs(ChoicePref_Correct(select_for_choicepref(:,st)&monkey2,j,st)), 'Tail', 'left');
            [~,p12]= ttest(abs(ChoicePref_Error(select_for_choicepref(:,st),j,st)), abs(ChoicePref_Correct(select_for_choicepref(:,st),j,st)), 'Tail', 'left');
            
        h= LinearCorrelation({
            ChoicePref_Error(monkey1&select_for_choicepref(:,st), j, st);
            ChoicePref_Error(monkey1&~select_for_choicepref(:,st),j, st);
            ChoicePref_Error(monkey2&select_for_choicepref(:,st),j, st);
            ChoicePref_Error(monkey2&~select_for_choicepref(:,st), j,st);
            },...
            {
            ChoicePref_Correct(monkey1&select_for_choicepref(:,st),j, st);
            ChoicePref_Correct(monkey1&~select_for_choicepref(:,st),j, st);
            ChoicePref_Correct(monkey2&select_for_choicepref(:,st), j,st);
            ChoicePref_Correct(monkey2&~select_for_choicepref(:,st), j,st);
            },...
            'Ylabel', 'Correct ChoicePref', 'Xlabel', 'Error ChoicePref',...
            'FaceColors', {condition_colors(st,:), [1 1 1], condition_colors(st,:), [1 1 1]},...
            'EdgeColors', {'k'}, 'Markers',{'o','o','^','^'},'MarkerSize', marker_size,...
            'LineStyles',{'k-'}, 'figN', figN, 'Axes',ha,...
            'SameScale',1, 'AxisSquare',1, 'Diagonal',1);
                
        delete([h.group.line]);
        
                % Show individual cell selected from the figure. HH20150424
        h_line = plot(ChoicePref_Error(:, j, st),ChoicePref_Correct(:,j,st),'visible','off'); hold on;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_bottom_line});

        
        text(0.1, -0.35-(st*0.2-0.1),...
            ['n= ',num2str(sum(select_for_choicepref(:,st)&monkey1)),' ', num2str(sum(select_for_choicepref(:,st)&monkey2)),' ', num2str(sum(select_for_choicepref(:,st)))],'color', condition_colors(st,:));
        text(0.1, -0.35-st*0.2, ['p: ' num2str(roundn(p1,-4)) '  ' num2str(roundn(p2,-4)) '  ' num2str(roundn(p12,-4))], 'color', condition_colors(st,:))

        axis([-1.1 1.1 -1.1 1.1]);
        xticks([-1 -0.5 0 0.5 1]); yticks([-1 -0.5 0 0.5 1]);
        
        legend off;
        figN = figN+1;
        
        SetFigure();
        end
    end


    function f3p2p5(debug)  % ChoiceDiv: Correct and Error
        if debug
            dbstack;
            keyboard;
        end
        j = 1;
        select_for_choicepref =  (squeeze(Choice_pref_p_value_all(:,:,j)) < 0.01)';
        
        monkeys = xls_num{1}(:,header.Monkey);  % HH20160613
        monkey1 = monkeys == 15; monkey1 = monkey1(select_bottom_line);
        monkey2 = monkeys == 13; monkey2 = monkey2(select_bottom_line);
        
        set(figure(figN),'name','Selected 300ms of ChoiceDiv Correct vs. Error'); clf; hold on; axis off;
        ha = axes('position', [0.1 0.1 0.8 0.8]);
        
        for st = 1:3
            
            % paired t-test
            [~,p1] = ttest(abs(mean(ChoiceDiv_Error{j}(select_for_choicepref(:,st)&monkey1,rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2),st),2)),...
                abs(mean(ChoiceDiv_Correct{j}(select_for_choicepref(:,st)&monkey1,rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2),st),2)), 'Tail', 'left');
            [~,p2] = ttest(abs(mean(ChoiceDiv_Error{j}(select_for_choicepref(:,st)&monkey2,rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2),st),2)),...
                abs(mean(ChoiceDiv_Correct{j}(select_for_choicepref(:,st)&monkey2,rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2),st),2)), 'Tail', 'left');
            [~,p12] = ttest(abs(mean(ChoiceDiv_Error{j}(select_for_choicepref(:,st),rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2),st),2)),...
                abs(mean(ChoiceDiv_Correct{j}(select_for_choicepref(:,st),rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2),st),2)), 'Tail', 'left');
            
            h= LinearCorrelation({
                mean(ChoiceDiv_Error{j}(monkey1&select_for_choicepref(:,st), rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2), st),2);
                mean(ChoiceDiv_Error{j}(monkey2&select_for_choicepref(:,st), rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2), st),2)                %             ChoicePref_Error(monkey2&~select_for_choicepref(:,st),j, st);
                },...
                {
                mean(ChoiceDiv_Correct{j}(monkey1&select_for_choicepref(:,st), rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2), st),2);
                mean(ChoiceDiv_Correct{j}(monkey2&select_for_choicepref(:,st), rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2), st),2)                %             ChoicePref_Error(monkey2&~select_for_choicepref(:,st),j, st);
                },...
                'Ylabel', 'Correct ChoiceDiv', 'Xlabel', 'Error ChoiceDiv',...
                'FaceColors', {condition_colors(st,:), condition_colors(st,:)},...
                'EdgeColors', {'k'}, 'Markers',{'o','^'},'MarkerSize', marker_size,...
                'LineStyles',{'k-'}, 'figN', figN, 'Axes',ha,...
                'SameScale',1, 'AxisSquare',1, 'Diagonal',1);
            
            delete([h.group.line]);
            
            text(0.1, -0.1-(st*0.1-0.05),...
                ['n= ',num2str(sum(select_for_choicepref(:,st)&monkey1)),' ', num2str(sum(select_for_choicepref(:,st)&monkey2)),' ', num2str(sum(select_for_choicepref(:,st)))],'color', condition_colors(st,:));
            text(0.1, -0.1-st*0.1, ['p: ' num2str(roundn(p1,-4)) '  ' num2str(roundn(p2,-4)) '  ' num2str(roundn(p12,-4))], 'color', condition_colors(st,:))
        end
        axis([-0.6 1.1 -0.6 1.1]);
        legend off;
        figN = figN+1;
        
        SetFigure();

    end

    function f3p2p6(debug)  % Choice Preference: Difficult and Easy
        if debug
            dbstack;
            keyboard;
        end
        
        j = 1; 
        select_for_choicepref =  (squeeze(Choice_pref_p_value_all(:,:,j)) < 0.01)';
        
        monkeys = xls_num{1}(:,header.Monkey);  % HH20160613
        monkey1 = monkeys == 15; monkey1 = monkey1(select_bottom_line);
        monkey2 = monkeys == 13; monkey2 = monkey2(select_bottom_line);

        set(figure(figN),'name','ChoicePref Easy vs. Difficult'); clf; hold on; axis off;
        ha = axes('position', [0.1 0.1 0.8 0.8]);
        
        for st = 1:3
            
            % paired t-test
            [~,p1] = ttest(abs(ChoicePref_Difficult(select_for_choicepref(:,st)&monkey1,j,st)), abs(ChoicePref_Easy(select_for_choicepref(:,st)&monkey1,j,st)), 'Tail', 'left');
            [~,p2] = ttest(abs(ChoicePref_Difficult(select_for_choicepref(:,st)&monkey2,j,st)), abs(ChoicePref_Easy(select_for_choicepref(:,st)&monkey2,j,st)), 'Tail', 'left');
            [~,p12]= ttest(abs(ChoicePref_Difficult(select_for_choicepref(:,st),j,st)), abs(ChoicePref_Easy(select_for_choicepref(:,st),j,st)), 'Tail', 'left');
            
            h= LinearCorrelation({
                ChoicePref_Difficult(monkey1&select_for_choicepref(:,st), j, st);
                %             ChoicePref_Error(monkey1&~select_for_choicepref(:,st),j, st);
                ChoicePref_Difficult(monkey2&select_for_choicepref(:,st),j, st);
                %             ChoicePref_Error(monkey2&~select_for_choicepref(:,st),j, st);
                },...
                {
                ChoicePref_Easy(monkey1&select_for_choicepref(:,st),j, st);
                %             ChoicePref_Correct(monkey1&~select_for_choicepref(:,st),j, st);
                ChoicePref_Easy(monkey2&select_for_choicepref(:,st), j,st);
                %             ChoicePref_Correct(monkey2&~select_for_choicepref(:,st), j,st);
                },...
                'Ylabel', 'Easy ChoicePref', 'Xlabel', 'Difficult ChoicePref',...
                'FaceColors', {condition_colors(st,:), condition_colors(st,:)},...
                'EdgeColors', {'k'}, 'Markers',{'o','^'},'MarkerSize', marker_size,...
                'LineStyles',{'k-'}, 'figN', figN, 'Axes',ha,...
                'SameScale',1, 'AxisSquare',1, 'Diagonal',1);
           
            delete([h.group.line]);
            
            text(0.1, -0.2-(st*0.2-0.1),...
                ['n= ',num2str(sum(select_for_choicepref(:,st)&monkey1)),' ', num2str(sum(select_for_choicepref(:,st)&monkey2)),' ', num2str(sum(select_for_choicepref(:,st)))],'color', condition_colors(st,:));
            text(0.1, -0.2-st*0.2, ['p: ' num2str(roundn(p1,-4)) '  ' num2str(roundn(p2,-4)) '  ' num2str(roundn(p12,-4))], 'color', condition_colors(st,:))
        end
        axis([-1.1 1.1 -1.1 1.1]); 
        xticks([-1 -0.5 0 0.5 1]); yticks([-1 -0.5 0 0.5 1])
        legend off;
        figN = figN+1;
        
        SetFigure();
    end


    function f3p2p7(debug)  % ChoiceDiv: Difficult and Easy
        if debug
            dbstack;
            keyboard;
        end
        
        j = 1; 
        select_for_choicepref =  (squeeze(Choice_pref_p_value_all(:,:,j)) < 0.01)';
        
        monkeys = xls_num{1}(:,header.Monkey);  % HH20160613
        monkey1 = monkeys == 15; monkey1 = monkey1(select_bottom_line);
        monkey2 = monkeys == 13; monkey2 = monkey2(select_bottom_line);

        set(figure(figN),'name','Selected 300ms of ChoiceDiv Easy vs. Difficult'); clf; hold on; axis off;
        ha = axes('position', [0.1 0.1 0.8 0.8]);
        
        for st = 1:3
            
            % paired t-test
            [~,p1] = ttest(abs(mean(ChoiceDiv_Difficult{j}(select_for_choicepref(:,st)&monkey1,rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2),st),2)),...
                abs(mean(ChoiceDiv_Easy{j}(select_for_choicepref(:,st)&monkey1,rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2),st),2)), 'Tail', 'left');
            [~,p2] = ttest(abs(mean(ChoiceDiv_Difficult{j}(select_for_choicepref(:,st)&monkey2,rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2),st),2)),...
                abs(mean(ChoiceDiv_Easy{j}(select_for_choicepref(:,st)&monkey2,rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2),st),2)), 'Tail', 'left');
            [~,p12] = ttest(abs(mean(ChoiceDiv_Difficult{j}(select_for_choicepref(:,st),rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2),st),2)),...
                abs(mean(ChoiceDiv_Easy{j}(select_for_choicepref(:,st),rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2),st),2)), 'Tail', 'left');
            
            format shortE;

        h= LinearCorrelation({
            mean(ChoiceDiv_Difficult{j}(monkey1&select_for_choicepref(:,st), rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2), st),2);
            mean(ChoiceDiv_Difficult{j}(monkey2&select_for_choicepref(:,st), rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2), st),2)                %             ChoicePref_Error(monkey2&~select_for_choicepref(:,st),j, st);
            },...
            {
            mean(ChoiceDiv_Easy{j}(monkey1&select_for_choicepref(:,st), rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2), st),2);
            mean(ChoiceDiv_Easy{j}(monkey2&select_for_choicepref(:,st), rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2), st),2)                %             ChoicePref_Error(monkey2&~select_for_choicepref(:,st),j, st);
            },...
            'Ylabel', 'Easy ChoiceDiv', 'Xlabel', 'Difficult ChoiceDiv',...
            'FaceColors', {condition_colors(st,:), condition_colors(st,:)},...
            'EdgeColors', {'k'}, 'Markers',{'o','^'},'MarkerSize', marker_size,...
            'LineStyles',{'k-'}, 'figN', figN, 'Axes',ha,...
            'SameScale',1, 'AxisSquare',1, 'Diagonal',1);
        
            delete([h.group.line]);
            
                        % Show individual cell selected from the figure. HH20150424
            h_line = plot(mean(ChoiceDiv_Difficult{j}(select_for_choicepref(:,st), rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2), st),2),...
            mean(ChoiceDiv_Easy{j}(select_for_choicepref(:,st), rate_ts{j}>tWin4CD{st,1}(1) & rate_ts{j}<tWin4CD{st,1}(2), st),2),'visible','off'); hold on;
            set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_bottom_line});

            
            text(0.1, -0.1-(st*0.1-0.05),...
                ['n= ',num2str(sum(select_for_choicepref(:,st)&monkey1)),' ', num2str(sum(select_for_choicepref(:,st)&monkey2)),' ', num2str(sum(select_for_choicepref(:,st)))],'color', condition_colors(st,:));
            text(0.1, -0.1-st*0.1, ['p: ' num2str(p1) '  ' num2str(p2) '  ' num2str(p12)], 'color', condition_colors(st,:))
%             text(0.1, -0.1-st*0.1, ['p: ' num2str(roundn(p1,-4)) '  ' num2str(roundn(p2,-4)) '  ' num2str(roundn(p12,-4))], 'color', condition_colors(st,:))
        end
        axis([-0.6 1.1 -0.6 1.1]);
        legend off;
        figN = figN+1;
        
        SetFigure();

    end


    function f3p1(debug)      % Correlations 1. Mem-sac v.s. CDiv
        if debug
            dbstack;
            keyboard;
        end
        
        
        monkeys = xls_num{1}(:,header.Monkey);  % HH20160613
        
        % ------- 1. Stimulus period vs. MemSac --------
        j = 1;
        set(figure(figN), 'name', 'Selected 300ms CDiv with MemSac'); clf; hold on; figN=figN+1;
        for c = 1:3 % condition
            for mon = monkey_included_for_analysis
                m = find(monkey_included_for_analysis == mon);
                
                ha = plot(MemSac_indicator(monkeys==mon & select_bottom_line),...
                    abs(mean(ChoiceDiv_All{j}(monkeys==mon & select_bottom_line, rate_ts{j} > tWin4CD{c,1}(1) & rate_ts{j} < tWin4CD{c,1}(2),c),2)),...
                    monkey_marker{m}, 'Color',condition_colors(c,:), 'MarkerSize',8, 'LineWidth',0.75);
                
                % Orthogonal regression for line
                fitPerp = regress_perp(MemSac_indicator(monkeys==mon & select_bottom_line),...
                    abs(mean(ChoiceDiv_All{j}(monkeys==mon & select_bottom_line, rate_ts{j} > tWin4CD{c,1}(1) & rate_ts{j} < tWin4CD{c,1}(2),c),2)),...
                    1, 0.05);
                linPara(1) = fitPerp.k;
                linPara(2) = fitPerp.b;
                linParaSE = fitPerp.kInterval;
                
                xx = linspace(min(MemSac_indicator(monkeys==mon & select_bottom_line)),...
                    max(MemSac_indicator(monkeys==mon & select_bottom_line)), 1000);
                %                 yy = polyval(linPara, xx);
                yy = linPara(1) * xx + linPara(2);
                
                axes = gca;
                xx = xx(yy > axes.YLim(1) & yy< axes.YLim(2));
                yy = yy(yy > axes.YLim(1) & yy< axes.YLim(2));
                
                plot(xx,yy,monkey_line{m}, 'color', condition_colors(c,:), 'linew', 2);
                
                % Pearson correlation
                [r, p] = corr(MemSac_indicator(monkeys==mon & select_bottom_line),...
                    abs(mean(ChoiceDiv_All{j}(monkeys==mon & select_bottom_line, rate_ts{j} > tWin4CD{c,1}(1) & rate_ts{j} < tWin4CD{c,1}(2),c),2)),...
                    'type', 'Spearman', 'row','complete');
                text(0.5, 1-0.03*(c+(m-1)*3),...
                    ['r^{2}= ' num2str(r^2) ' p: ' num2str(p) ' k= ' num2str(linPara(1)) ' ( ' num2str(linParaSE) ' )'],...
                    'color', condition_colors(c,:))
                
                % Annotate tcells
                plot(MemSac_indicator(monkeys==mon & select_tcells),...
                    abs(mean(ChoiceDiv_All{j}(monkeys==mon & select_tcells,rate_ts{j} > tWin4CD{c,1}(1) & rate_ts{j} < tWin4CD{c,1}(2),c),2)),...
                    monkey_marker{m}, 'Color',condition_colors(c,:), 'MarkerSize',8, 'LineWidth',0.75, 'MarkerFaceColor',condition_colors(c,:));
            end
        end
        xlabel('MemSac (DDI)'); ylabel('abs(Stimuli CDiv) ');
        SetFigure();
        
        
        % ---------- 2. Pre-saccade vs. Memsac ------------
        j = 2;
        set(figure(figN), 'name', 'Pre-sac CDiv with MemSac'); clf; hold on; figN=figN+1;
        for c = 1:3 % condition
            for mon = monkey_included_for_analysis
                m = find(monkey_included_for_analysis == mon);
                
                ha = plot(MemSac_indicator(monkeys==mon & select_bottom_line),...
                    abs(mean(ChoiceDiv_All{j}(monkeys==mon & select_bottom_line, rate_ts{j} > -300 & rate_ts{j} < 0,c),2)),...
                    monkey_marker{m}, 'Color',condition_colors(c,:), 'MarkerSize',8, 'LineWidth',0.75);
                
                % Orthogonal regression for line
                fitPerp = regress_perp(MemSac_indicator(monkeys==mon & select_bottom_line),...
                    abs(mean(ChoiceDiv_All{j}(monkeys==mon & select_bottom_line, rate_ts{j} > -300 & rate_ts{j} < 0,c),2)),...
                    1, 0.05);
                linPara(1) = fitPerp.k;
                linPara(2) = fitPerp.b;
                linParaSE = fitPerp.kInterval;
                
                xx = linspace(min(MemSac_indicator(monkeys==mon & select_bottom_line)),...
                    max(MemSac_indicator(monkeys==mon & select_bottom_line)), 1000);
                %                 yy = polyval(linPara, xx);
                yy = linPara(1) * xx + linPara(2);
                
                axes = gca;
                xx = xx(yy > axes.YLim(1) & yy< axes.YLim(2));
                yy = yy(yy > axes.YLim(1) & yy< axes.YLim(2));
                
                plot(xx,yy,monkey_line{m}, 'color', condition_colors(c,:), 'linew', 2);
                
                % Pearson correlation
                [r, p] = corr(MemSac_indicator(monkeys==mon & select_bottom_line),...
                    abs(mean(ChoiceDiv_All{j}(monkeys==mon & select_bottom_line, rate_ts{j} > -300 & rate_ts{j} < 0,c),2)),...
                    'type', 'Spearman', 'row','complete');
                text(0.5, 1-0.03*(c+(m-1)*3),...
                    ['r^{2}= ' num2str(r^2) ' p: ' num2str(p) ' k= ' num2str(linPara(1)) ' ( ' num2str(linParaSE) ' )'],...
                    'color', condition_colors(c,:))
                
                % Annotate tcells
                plot(MemSac_indicator(monkeys==mon & select_tcells),...
                    abs(mean(ChoiceDiv_All{j}(monkeys==mon & select_tcells, rate_ts{j} > -300 & rate_ts{j} < 0,c),2)),...
                    monkey_marker{m}, 'Color',condition_colors(c,:), 'MarkerSize',8, 'LineWidth',0.75, 'MarkerFaceColor',condition_colors(c,:));
            end
        end
        xlabel('MemSac (DDI)'); ylabel('abs(Pre-sac CDiv)');
        SetFigure();
        
        
        % ---------- 3. Feedback vs. Memsac ------------
        j = 3;
        set(figure(figN), 'name', 'Feedback CDiv with MemSac'); clf; hold on; figN=figN+1;
        for c = 1:3 % condition
            for mon = monkey_included_for_analysis
                m = find(monkey_included_for_analysis == mon);
                
                ha = plot(MemSac_indicator(monkeys==mon & select_bottom_line),...
                    abs(mean(ChoiceDiv_All{j}(monkeys==mon & select_bottom_line, rate_ts{j} > 50 & rate_ts{j} < 350,c),2)),...
                    monkey_marker{m}, 'Color',condition_colors(c,:), 'MarkerSize',8, 'LineWidth',0.75);
                
                % Orthogonal regression for line
                fitPerp = regress_perp(MemSac_indicator(monkeys==mon & select_bottom_line),...
                    abs(mean(ChoiceDiv_All{j}(monkeys==mon & select_bottom_line, rate_ts{j} > 50 & rate_ts{j} < 350,c),2)),...
                    1, 0.05);
                linPara(1) = fitPerp.k;
                linPara(2) = fitPerp.b;
                linParaSE = fitPerp.kInterval;
                
                xx = linspace(min(MemSac_indicator(monkeys==mon & select_bottom_line)),...
                    max(MemSac_indicator(monkeys==mon & select_bottom_line)), 1000);
                %                 yy = polyval(linPara, xx);
                yy = linPara(1) * xx + linPara(2);
                
                axes = gca;
                xx = xx(yy > axes.YLim(1) & yy< axes.YLim(2));
                yy = yy(yy > axes.YLim(1) & yy< axes.YLim(2));
                
                plot(xx,yy,monkey_line{m}, 'color', condition_colors(c,:), 'linew', 2);
                
                % Pearson correlation
                [r, p] = corr(MemSac_indicator(monkeys==mon & select_bottom_line),...
                    abs(mean(ChoiceDiv_All{j}(monkeys==mon & select_bottom_line, rate_ts{j} > 50 & rate_ts{j} < 350,c),2)),...
                    'type', 'Spearman', 'row','complete');
                text(0.5, 1-0.03*(c+(m-1)*3),...
                    ['r^{2}= ' num2str(r^2) ' p: ' num2str(p) ' k= ' num2str(linPara(1)) ' ( ' num2str(linParaSE) ' )'],...
                    'color', condition_colors(c,:))
                
                % Annotate tcells
                plot(MemSac_indicator(monkeys==mon & select_tcells),...
                    abs(mean(ChoiceDiv_All{j}(monkeys==mon & select_tcells, rate_ts{j} > 50 & rate_ts{j} < 350,c),2)),...
                    monkey_marker{m}, 'Color',condition_colors(c,:), 'MarkerSize',8, 'LineWidth',0.75, 'MarkerFaceColor',condition_colors(c,:));
            end
        end
        xlabel('MemSac (DDI)'); ylabel('abs(Feedback CDiv)');
        SetFigure();
        
        
        nHist = 10;
        sig_marker_size = 10;
        % =================== Correlations 2. Mem-sac v.s. CP
        
        j = 1; % Not much time-sensitive, so I use j = 2 here.
        %         CP_interest = [-300 0];  % min and max interval. HH20180609
        CP_interest = tWin4CD;   % Changed byZZ @20230628
        
        %  CP_interval = CP_ts{j} == CP_interest ;
        
        for k = 1:3
            CP_interval{k} = find(abs(CP_ts{j} - CP_interest{k,1}(1)) == min(abs(CP_ts{j} - CP_interest{k,1}(1))),1) : ...
                find(abs(CP_ts{j} - CP_interest{k,1}(2)) == min(abs(CP_ts{j} - CP_interest{k,1}(2))));
            CP_NS{k} = find_bottom_line(any(CP_p{j}(select_bottom_line,CP_interval{k},k) > 0.05,2));
            CP_S{k} = find_bottom_line(any(CP_p{j}(select_bottom_line,CP_interval{k},k) <= 0.05,2));
        end
        
        
        h = LinearCorrelation(...
            {   MemSac_indicator(CP_NS{1}),...
            MemSac_indicator(CP_S{1}),...
            MemSac_indicator(CP_NS{2}),...
            MemSac_indicator(CP_S{2}),...
            MemSac_indicator(CP_NS{3}),...
            MemSac_indicator(CP_S{3}),...
            },...
            {  mean(CP{j}(CP_NS{1},CP_interval{1},1),2),...
            mean(CP{j}(CP_S{1},CP_interval{1},1),2),...
            mean(CP{j}(CP_NS{2},CP_interval{2},2),2),...
            mean(CP{j}(CP_S{2},CP_interval{2},2),2),...
            mean(CP{j}(CP_NS{3},CP_interval{3},3),2),...
            mean(CP{j}(CP_S{3},CP_interval{3},3),2)},...
            'CombinedIndex',[3 12 48],...
            'MethodOfCorr','Spearman', 'FittingMethod', 2,...
            'Xlabel',MemSac_indicator_txt,'Ylabel','CP at selected stimulus 300ms ',...
            'FaceColors',{'none',condition_colors(1,:),'none',condition_colors(2,:),'none',condition_colors(3,:)},...
            'EdgeColors',{condition_colors(1,:),'none',condition_colors(2,:),'none',condition_colors(3,:),'none'}, 'Markers',{'o'},...
            'LineStyles',{'b:','b-','r:','r-','g:','g-','b-','r-','g-'},'MarkerSize',marker_size,...
            'figN',figN,'XHist',nHist,'YHist',nHist); figN = figN + 1;
        set(gcf,'name',['j = ' num2str(j), 'MemSac vs. CP']);
        delete([h.group(1:6).line]);
        plot(xlim,[0.5 0.5],'k--');         SetFigure(20);
        
        
        % Annotate tcells
        for k = 1:3
            plot(MemSac_indicator(select_tcells),mean(CP{j}(select_tcells,CP_interval{k},k),2),...
                '+','markersize',tcell_cross_size,'color','k','linew',2);
        end
        
        % ========== CD and CP ========
        h = LinearCorrelation({  mean(ChoiceDiv_All{j}(CP_NS{1},rate_ts{j} >tWin4CD{1,1}(1) & rate_ts{j} < tWin4CD{1,1}(2),1),2) / 2 + 0.5,...
            mean(ChoiceDiv_All{j}(CP_S{1},rate_ts{j} > tWin4CD{1,1}(1) & rate_ts{j} < tWin4CD{1,1}(2),1),2) / 2 + 0.5,...
            mean(ChoiceDiv_All{j}(CP_NS{2},rate_ts{j} > tWin4CD{2,1}(1) & rate_ts{j} < tWin4CD{2,1}(2),2),2)/ 2 + 0.5,...
            mean(ChoiceDiv_All{j}(CP_S{2},rate_ts{j} > tWin4CD{2,1}(1) & rate_ts{j} < tWin4CD{2,1}(2),2),2)/ 2 + 0.5,...
            mean(ChoiceDiv_All{j}(CP_NS{3},rate_ts{j} > tWin4CD{3,1}(1) & rate_ts{j} < tWin4CD{3,1}(2),3),2)/ 2 + 0.5,...
            mean(ChoiceDiv_All{j}(CP_S{3},rate_ts{j} > tWin4CD{3,1}(1) & rate_ts{j} < tWin4CD{3,1}(2),3),2)/ 2 + 0.5},...
            {  mean(CP{j}(CP_NS{1},CP_interval{1},1),2),...
            mean(CP{j}(CP_S{1},CP_interval{1},1),2),...
            mean(CP{j}(CP_NS{2},CP_interval{2},2),2),...
            mean(CP{j}(CP_S{2},CP_interval{2},2),2),...
            mean(CP{j}(CP_NS{3},CP_interval{3},3),2),...
            mean(CP{j}(CP_S{3},CP_interval{3},3),2)},...
            'CombinedIndex',[3 12 48],...
            'MethodOfCorr','Spearman', 'FittingMethod', 2,...
            'Xlabel','Choice Div (pre)/2 + 0.5','Ylabel','CP at selected stimulus 300ms ',...
            'FaceColors',{'none',condition_colors(1,:),'none',condition_colors(2,:),'none',condition_colors(3,:)},...
            'EdgeColors',{condition_colors(1,:),'none',condition_colors(2,:),'none',condition_colors(3,:),'none'}, 'Markers',{'o'},...
            'LineStyles',{'b:','b-','r:','r-','g:','g-','b-','r-','g-'},'MarkerSize',marker_size,...
            'figN',figN,'XHist',nHist,'YHist',nHist,'SameScale',1); figN = figN + 1;
        delete([h.group(1:6).line]);
        plot(xlim,[0.5 0.5],'k--'); plot([0.5 0.5],ylim,'k--');         SetFigure(20);
        
        
        set(gcf,'name',['j = ' num2str(j) 'CD vs. CP']);
        % Annotate tcells
        for k = 1:3
            plot(.5 + .5* mean(ChoiceDiv_All{j}(select_tcells,rate_ts{j} > tWin4CD{k,1}(1) & rate_ts{j} < tWin4CD{k,1}(2),k),2),...
                mean(CP{j}(select_tcells,CP_interval{k},k),2),...
                '+','markersize',tcell_cross_size,'color','k','linew',2);
        end
    end

    function f3p1p2(debug)      % Mem-sac v.s. CPref
        if debug
            dbstack;
            keyboard;
        end
        % ------------ Memsac DDI and Choice Preference (single modality)
        
        monkeys = xls_num{1}(:,header.Monkey);  % HH20160613
        monkey1 = monkeys == 15; monkey1 = monkey1(select_bottom_line)';
        monkey2 = monkeys == 13; monkey2 = monkey2(select_bottom_line)';
        
        memsac_sig = MemSac_indicator_p(select_bottom_line)' < 0.05;
        
        for k = 1:3
            
            tt = 1; % From stim on to stim off
            
            cpref_sig = Choice_pref_p_value_all(k,:,tt) < 0.01;
            
            MemSac_DDI_selected = MemSac_indicator(select_bottom_line);
            
            h = LinearCorrelation({
                MemSac_DDI_selected(monkey1 & ~or(cpref_sig, memsac_sig));
                MemSac_DDI_selected(monkey1 & xor(cpref_sig, memsac_sig));
                MemSac_DDI_selected(monkey1 & and(cpref_sig, memsac_sig));
                MemSac_DDI_selected(monkey2 & ~or(cpref_sig, memsac_sig));
                MemSac_DDI_selected(monkey2 & xor(cpref_sig, memsac_sig));
                MemSac_DDI_selected(monkey2 & and(cpref_sig, memsac_sig));
                },...
                {
                abs(Choice_pref_all(k,monkey1 & ~or(cpref_sig, memsac_sig),tt)) ;
                abs(Choice_pref_all(k,monkey1 & xor(cpref_sig, memsac_sig),tt));
                abs(Choice_pref_all(k,monkey1 & and(cpref_sig, memsac_sig),tt));
                abs(Choice_pref_all(k,monkey2 & ~or(cpref_sig, memsac_sig),tt)) ;
                abs(Choice_pref_all(k,monkey2 & xor(cpref_sig, memsac_sig),tt)) ;
                abs(Choice_pref_all(k,monkey2 & and(cpref_sig, memsac_sig),tt));
                },...
                'CombinedIndex',[7 56 54 63],...
                'Xlabel',MemSac_indicator_txt,'Ylabel','abs(Choice preference)',...
                'FaceColors',{'none',condition_colors(k,:)*0.2 + [0.8 0.8 0.8],condition_colors(k,:)},'Markers',{'o','o','o','^','^','^'},...
                'EdgeColors',mat2cell(repmat(condition_colors(k,:),6,1),ones(6,1)),...
                'LineStyles',{'k-','k:','k--','k-.'},'MarkerSize',15,...
                'figN',figN,'XHist',20,'YHist',20,...
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,...
                'Method','Spearman','FittingMethod',2); figN = figN + 1;
            
            
            % Show individual cell selected from the figure. HH20150424
            h_line = plot(MemSac_DDI_selected,Choice_pref_all(k,:,tt),'visible','off'); hold on;
            set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_bottom_line});
            
            
            ylim([0 1]);
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
            %             delete([h.group(1:4).line]);
            % Annotate tcells
%             plot(MemSac_DDI_selected(select_tcells),abs(Choice_pref_all(k,select_tcells,tt)),...
%                 '+','markersize',tcell_cross_size,'color','k','linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            h_line = plot(MemSac_DDI_selected(:),abs(Choice_pref_all(k,:,tt)),'visible','off'); hold on;
            set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_bottom_line});
            
            % Override histogram
            cla(h.ax_xhist);
            HistComparison( {MemSac_DDI_selected(memsac_sig), MemSac_DDI_selected(~memsac_sig)},...
                'EdgeColors',{'k','k'}, 'FaceColors',{'k', 'none'},...
                'MeanType','Median','Axes',h.ax_xhist,'XCenters',0:0.05:1,'TTest',0);
            h.ax_xhist.XLim = xlim(h.ax_raw);
            
            cla(h.ax_yhist);
            HistComparison( {abs(Choice_pref_all(k,cpref_sig,tt))',abs(Choice_pref_all(k, ~cpref_sig,tt))'},...
                'EdgeColors',{'k','k'}, 'FaceColors',{condition_colors(k,:), 'none'},...
                'MeanType','Median','Axes',h.ax_yhist,'XCenters',0:0.05:1,'TTest',0);
            h.ax_yhist.XLim = ylim(h.ax_raw);
            
        end
    end
    function f3p2(debug)      % Correlations 2. Choice preference v.s. modality preference (Anne Fig.2f-h; Fig.3a)
        if debug
            dbstack;
            keyboard;
        end
        
        %% ---------------- 1 Choice pref v.s. modality pref (Anne Fig.3a)-------------------
        %         tt = 3;
        tt = 1;
        
        monkeys = xls_num{1}(:,header.Monkey);  % HH20160613
        monkey1 = monkeys == 15; monkey1 = monkey1(select_bottom_line)';
        monkey2 = monkeys == 13; monkey2 = monkey2(select_bottom_line)';
        
        for k = 1:3
            set(figure(figN),'name','CDiv vs. MDiv','pos',[17 514 1151 449]);
            cpref_sig = Choice_pref_p_value_all(k,:,tt) < 0.01;
            mpref_sig = Modality_pref_p_value_all(1,:,tt) < 0.01;
            
            h = LinearCorrelation({
                (Modality_pref_all(1, monkey1 & ~cpref_sig & ~mpref_sig,tt)) ;
                (Modality_pref_all(1, monkey2 & ~cpref_sig & ~mpref_sig,tt)) ;
                (Modality_pref_all(1, monkey1 & xor(cpref_sig , mpref_sig),tt));
                (Modality_pref_all(1, monkey2 & xor(cpref_sig , mpref_sig),tt));
                (Modality_pref_all(1, monkey1 & cpref_sig & mpref_sig,tt));...
                (Modality_pref_all(1, monkey2 & cpref_sig & mpref_sig,tt))
                },...
                {
                (Choice_pref_all(k,monkey1 & ~cpref_sig & ~mpref_sig,tt)) ;
                (Choice_pref_all(k,monkey2 & ~cpref_sig & ~mpref_sig,tt)) ;
                (Choice_pref_all(k,monkey1 & xor(cpref_sig , mpref_sig),tt)) ;
                (Choice_pref_all(k,monkey2 & xor(cpref_sig , mpref_sig),tt)) ;
                (Choice_pref_all(k,monkey1 & cpref_sig & mpref_sig,tt)) ;...
                (Choice_pref_all(k,monkey2 & cpref_sig & mpref_sig,tt))
                },...
                'CombinedIndex',[21 42 63],...
                'Ylabel','Choice preference (pre)','Xlabel','Modality preference (pre)',...
                'FaceColors',{'w','w',condition_colors(k,:)*0.2 + [0.8 0.8 0.8],condition_colors(k,:)*0.2 + [0.8 0.8 0.8],condition_colors(k,:),condition_colors(k,:)},'Markers',{'o','^'},...
                'EdgeColors',mat2cell(repmat(condition_colors(k,:),6,1),ones(6,1)),...
                'LineStyles',{'k:','k--','k-'},'MarkerSize',marker_size,...
                'figN',figN,'XHist',20,'YHist',20,...
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',1,...
                'Method','Spearman','FittingMethod',2); figN = figN + 1;
            
            %             delete([h.group(1:6).line h.diag]);
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');       SetFigure(15);
            axis([-1 1 -1 1])
            set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1);
            
            
            % Annotate tcells
            plot(Modality_pref_all(1,select_tcells(select_bottom_line),tt),...
                Choice_pref_all(k,select_tcells(select_bottom_line),tt),'+','markersize',tcell_cross_size,'color','k','linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            h_line = plot(Modality_pref_all(1,:,tt),Choice_pref_all(k,:,tt),'visible','off'); hold on;
            set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_cpref_mpref});
            
            % Override histogram
            cla(h.ax_xhist);
            HistComparison( {Modality_pref_all(1,mpref_sig,tt)', Modality_pref_all(1,~mpref_sig,tt)'},...
                'EdgeColors',{'k','k'}, 'FaceColors',{'k', 'none'},...
                'MeanType','Mean','Axes',h.ax_xhist,'XCenters',-1:0.1:1,'TTest',0);
            h.ax_xhist.XLim = xlim(h.ax_raw);
            
            cla(h.ax_yhist);
            HistComparison( {Choice_pref_all(k,cpref_sig,tt)', Choice_pref_all(k,~cpref_sig,tt)'},...
                'EdgeColors',{'k','k'}, 'FaceColors',{'k', 'none'},...
                'MeanType','Mean','Axes',h.ax_yhist,'XCenters',-1:0.1:1,'TTest',0);
            h.ax_yhist.XLim = ylim(h.ax_raw);
            
        end
    end

    function f3p2p2(debug)      % Correlations 2. Choice preference between modalities
        if debug
            dbstack;
            keyboard;
        end
        %% ---------------- 2 Choice pref (vest and visual) -------------------
        
        %         tt = 3;
        tt = 1;
        
        monkeys = xls_num{1}(:,header.Monkey);  % HH20160613
        monkey1 = monkeys == 15; monkey1 = monkey1(select_bottom_line)';
        monkey2 = monkeys == 13; monkey2 = monkey2(select_bottom_line)';
        
        
        set(figure(figN),'name','CDiv (visual) vs. CDiv (vest)','pos',[17 514 1151 449]);clf; 
        
        cpref_sig_1 = Choice_pref_p_value_all(1,:,tt) <= 0.01;
        cpref_sig_2 = Choice_pref_p_value_all(2,:,tt) <= 0.01;
        
        % Abs() or not?
        Choice_pref_all_temp = (Choice_pref_all);
        
        h = LinearCorrelation({
            (Choice_pref_all_temp(2, monkey1 & ~cpref_sig_1 & ~cpref_sig_2,tt)) ;
            (Choice_pref_all_temp(2, monkey1 & xor(cpref_sig_1 , cpref_sig_2),tt));
            (Choice_pref_all_temp(2, monkey1 & cpref_sig_1 & cpref_sig_2,tt))
            (Choice_pref_all_temp(2, monkey2 & ~cpref_sig_1 & ~cpref_sig_2,tt)) ;
            (Choice_pref_all_temp(2, monkey2 & xor(cpref_sig_1 , cpref_sig_2),tt));
            (Choice_pref_all_temp(2, monkey2 & cpref_sig_1 & cpref_sig_2,tt))
            },...
            {
            (Choice_pref_all_temp(1,monkey1 & ~cpref_sig_1 & ~cpref_sig_2,tt)) ;
            (Choice_pref_all_temp(1,monkey1 & xor(cpref_sig_1 , cpref_sig_2),tt)) ;
            (Choice_pref_all_temp(1,monkey1 & cpref_sig_1 & cpref_sig_2,tt))
            (Choice_pref_all_temp(1,monkey2 & ~cpref_sig_1 & ~cpref_sig_2,tt)) ;
            (Choice_pref_all_temp(1,monkey2 & xor(cpref_sig_1 , cpref_sig_2),tt)) ;
            (Choice_pref_all_temp(1,monkey2 & cpref_sig_1 & cpref_sig_2,tt))
            },...
            'CombinedIndex',[63],'PlotCombinedOnly', 1, ...
            'Ylabel','Vestibular choice preference','Xlabel','Visual choice preference',...
            'FaceColors',{'none',[0.8 0.8 0.8],'k'},'Markers',{'o', '^'},...
            'EdgeColors',{'k'},...
            'LineStyles',{'k-'},'MarkerSize',marker_size,...
            'figN',figN,'XHist',20,'YHist',20,'Method','Pearson','FittingMethod',2, ...
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',1); figN = figN + 1;
        
        %         delete([h.group(1:6).line]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');         SetFigure(20);
        
        %         % Annotate tcells
        %         plot(Choice_pref_all_temp(2,select_tcells(select_bottom_line),tt),Choice_pref_all_temp(1,select_tcells(select_bottom_line),tt),...
        %             '+','markersize',tcell_cross_size,'color',colors(2,:),'linew',2);
        
        axis([-1 1 -1 1])
        set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot(Choice_pref_all(2,:,tt),Choice_pref_all(1,:,tt),'visible','off'); hold on;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_cpref_mpref});
        
        % Override histogram
        cla(h.ax_xhist);
        HistComparison( {Choice_pref_all_temp(1,cpref_sig_1,tt)', Choice_pref_all_temp(1,~cpref_sig_1,tt)'},...
            'EdgeColors',{'k','k'}, 'FaceColors',{'k', 'none'},...
            'MeanType','Median','Axes',h.ax_xhist,'XCenters',-1:0.1:1,'TTest',0);
        h.ax_xhist.XLim = xlim(h.ax_raw);
        
        cla(h.ax_yhist);
        HistComparison( {Choice_pref_all_temp(2,cpref_sig_2,tt)', Choice_pref_all_temp(2,~cpref_sig_2,tt)'},...
            'EdgeColors',{'k','k'}, 'FaceColors',{'k', 'none'},...
            'MeanType','Median','Axes',h.ax_yhist,'XCenters',-1:0.1:1,'TTest',0);
        h.ax_yhist.XLim = ylim(h.ax_raw);
        
        
        clear Choice_pref_all_temp;
        
        
        %% ---------------- 3 Choice pref (single vs comb) -------------------
        two_face_colors = fliplr({'none',[0.8 0.8 1],condition_colors(1,:);
            'none',[1 0.8 0.8],condition_colors(2,:)});
        
        for k = 1:2  % Plot it separately
            
            
            set(figure(figN),'name','CDiv (single) vs. CDiv (comb)','pos',[17 514 1151 449]);clf; 
            
            cpref_sig_k = Choice_pref_p_value_all(k,:,tt) <= 0.01;
            %             cpref_sig_2 = Choice_pref_p_value_all(2,:,tt) < 0.05;
            cpref_sig_3 = Choice_pref_p_value_all(3,:,tt) <= 0.01;
            
            % Abs() or not?
            Choice_pref_all_temp = (Choice_pref_all);
            
            
            h = LinearCorrelation({
                (Choice_pref_all_temp(k, monkey1 & cpref_sig_k & cpref_sig_3,tt));
                (Choice_pref_all_temp(k, monkey1 & xor(cpref_sig_k , cpref_sig_3),tt));
                (Choice_pref_all_temp(k, monkey1 & ~cpref_sig_k & ~cpref_sig_3,tt)) ;
                (Choice_pref_all_temp(k, monkey2 & cpref_sig_k & cpref_sig_3,tt));
                (Choice_pref_all_temp(k, monkey2 & xor(cpref_sig_k , cpref_sig_3),tt));
                (Choice_pref_all_temp(k, monkey2 & ~cpref_sig_k & ~cpref_sig_3,tt)) ;
                },...
                {
                (Choice_pref_all_temp(3,monkey1 & cpref_sig_k & cpref_sig_3,tt)) ;
                (Choice_pref_all_temp(3,monkey1 & xor(cpref_sig_k, cpref_sig_3),tt)) ;
                (Choice_pref_all_temp(3,monkey1 & ~cpref_sig_k & ~cpref_sig_3,tt)) ;
                (Choice_pref_all_temp(3,monkey2 & cpref_sig_k & cpref_sig_3,tt));
                (Choice_pref_all_temp(3,monkey2 & xor(cpref_sig_k, cpref_sig_3),tt));
                (Choice_pref_all_temp(3,monkey2 & ~cpref_sig_k & ~cpref_sig_3,tt)) ;
                },...
                'CombinedIndex',[63],'PlotCombinedOnly', 1, ...
                'Ylabel','Combined choice preference','Xlabel','Single choice preference',...
                'FaceColors',two_face_colors(k,:),'Markers',{'o','^'},'EdgeColors',{'k'},...
                'LineStyles',{'k-'},'MarkerSize',marker_size,...
                'figN',figN,'XHist',20,'YHist',20,...
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',1,...
                'Method','Pearson','FittingMethod',2); figN = figN + 1;
            
            % delete([h.group(1:6).line]);
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');          SetFigure(20);
            axis([-1 1 -1 1]); set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1); axis square;
            
            % Annotate tcells
            %             plot((Choice_pref_all_temp(k,select_tcells(select_bottom_line),tt)),(Choice_pref_all_temp(3,select_tcells(select_bottom_line),tt)),...
            %                 '+','markersize',tcell_cross_size,'color','k','linew',2);
            %             plot((Choice_pref_all_temp(2,select_tcells(select_bottom_line),tt)),(Choice_pref_all_temp(3,select_tcells(select_bottom_line),tt)),...
            %                 '+','markersize',tcell_cross_size,'color','k','linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            h_line = plot([Choice_pref_all_temp(k,:,tt)],[Choice_pref_all_temp(3,:,tt)],'visible','off'); hold on;
            set([gca h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, [select_cpref_mpref],1});
            
            % Plotting histograms manually because we have to deal with different significant standards here. HH20180613
            cla(h.ax_xhist);
            HistComparison( {Choice_pref_all_temp(k,cpref_sig_k,tt)', Choice_pref_all_temp(k,~cpref_sig_k,tt)'},...
                'EdgeColors',{'k','k'}, 'FaceColors',{condition_colors(k,:), 'none'},...
                'MeanType','Median','Axes',h.ax_xhist,'XCenters',-1:0.1:1,'TTest',0);
            xlim([-1 1])
            
            cla(h.ax_yhist);
            HistComparison( {Choice_pref_all_temp(3,cpref_sig_3,tt)', Choice_pref_all_temp(3,~cpref_sig_3,tt)'},...
                'EdgeColors',{'k','k'}, 'FaceColors',{condition_colors(3,:), 'none'},...
                'MeanType','Median','Axes',h.ax_yhist,'XCenters',-1:0.1:1,'TTest',0);
            xlim([-1 1])
            
            clear Choice_pref_all_temp;
            axis tight;
        end
        
        %% ===  (Comb - visual) VS (Comb - vest) %HH20160830 ===
        set(figure(figN),'name','Cdiv(Comb - visual) VS Cdiv(Comb - vest)','pos',[17 514 1151 449]);clf; 
        
        % Abs() or not?
        Choice_pref_all_temp_comb_minus_vest = abs(Choice_pref_all(3,:,:))-abs(Choice_pref_all(1,:,:));
        Choice_pref_all_temp_comb_minus_vis = abs(Choice_pref_all(3,:,:))-abs(Choice_pref_all(2,:,:));
        
        %         cellTypes = [group_result.Waveform_broad];
        
        h = LinearCorrelation({
            (Choice_pref_all_temp_comb_minus_vest(1,monkey1,tt)) ;
            (Choice_pref_all_temp_comb_minus_vest(1,monkey2,tt)) ;
            %             (Choice_pref_all_temp_comb_minus_vest(1,monkey1,tt)) ;
            %             (Choice_pref_all_temp_comb_minus_vest(1,monkey2 & cellTypes,tt)) ;
            },...
            {
            (Choice_pref_all_temp_comb_minus_vis(1,monkey1,tt)) ;
            (Choice_pref_all_temp_comb_minus_vis(1,monkey2 ,tt)) ;
            %             (Choice_pref_all_temp_comb_minus_vis(1,monkey1 & cellTypes,tt)) ;
            %             (Choice_pref_all_temp_comb_minus_vis(1,monkey2 & cellTypes,tt)) ;
            },...
            'CombinedIndex',[3],'PlotCombinedOnly', 0, ...
            'Xlabel','Combined - Vest (Choice preference)','Ylabel','Combined - Visual',...
            'FaceColors',{'k'},'Markers',{'o','^'},'EdgeColors',{'k','k'},...
            'LineStyles',{'k:','k--','k-'},'MarkerSize',marker_size,...
            'figN',figN,'XHist',20,'YHist',20,'XHistStyle','stacked','YHistStyle','stacked',...
            'SameScale',1,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        
        %         delete([h.group(1:2).line]);
        axis([-1 1 -1 1]); set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1); axis square;
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');         SetFigure(20);
        
        %         % Annotate tcells
        %         plot((Choice_pref_all_temp_comb_minus_vest(1,select_tcells(select_bottom_line),tt)),...
        %             (Choice_pref_all_temp_comb_minus_vis(1,select_tcells(select_bottom_line),tt)),...
        %             '+','markersize',tcell_cross_size,'color',colors(2,:),'linew',2);
        %         %             plot((Choice_pref_all_temp(2,select_tcells(select_bottom_line),tt)),(Choice_pref_all_temp(3,select_tcells(select_bottom_line),tt)),...
        %         %                 '+','markersize',tcell_cross_size,'color','k','linew',2);
        % annotated by ZZ 20210111
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot([Choice_pref_all_temp_comb_minus_vest(1,:,tt)],[Choice_pref_all_temp_comb_minus_vis(1,:,tt)],'visible','off'); hold on;
        set([gca h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, [select_cpref_mpref],1});
        
        cla(h.ax_xhist);
        HistComparison( {Choice_pref_all_temp_comb_minus_vest(1,:,tt)'},...
            'EdgeColors',{'k'}, 'FaceColors',{'k'},...
            'MeanType','Median','Axes',h.ax_xhist,'XCenters',-1:0.1:1,'TTest',0);
        xlim([-1 1])
        cla(h.ax_yhist);
        HistComparison( {Choice_pref_all_temp_comb_minus_vis(1,:,tt)'},...
            'EdgeColors',{'k'}, 'FaceColors',{'k'},...
            'MeanType','Median','Axes',h.ax_yhist,'XCenters',-1:0.1:1,'TTest',0);
        xlim([-1 1])
        
        
        clear Choice_pref_all_temp;
        
        %% == 3-D Plotting for Comparison
        
        Choice_pref_all_temp = squeeze(Choice_pref_all(:,:,tt));
        cpref_sig = squeeze(Choice_pref_p_value_all(:,:,tt)) < 0.01;
        
        all_sig = all(cpref_sig);
        any_sig = any(cpref_sig) - all_sig;
        
        set(figure(figN), 'name', 'Choice Preference from All Conditions'); figN=figN+1; clf; 
        hold on; 
        plot3(Choice_pref_all_temp(1,monkey1), Choice_pref_all_temp(2,monkey1), Choice_pref_all_temp(3,monkey1),...
            'ko','MarkerSize',14)
        plot3(Choice_pref_all_temp(1,monkey2), Choice_pref_all_temp(2,monkey2), Choice_pref_all_temp(3,monkey2),...
            'k^','MarkerSize',14)
        
        % Annotate cells
        % All cells
        plot3(Choice_pref_all_temp(1,monkey1&all_sig), Choice_pref_all_temp(2,monkey1&all_sig),...
            Choice_pref_all_temp(3,monkey1& all_sig),'ko','MarkerSize',14,'MarkerFaceColor','k')
        plot3(Choice_pref_all_temp(1,monkey2&all_sig), Choice_pref_all_temp(2,monkey2&all_sig),...
            Choice_pref_all_temp(3,monkey2& all_sig),'k^','MarkerSize',14,'MarkerFaceColor','k')
        % Any cells
        plot3(Choice_pref_all_temp(1,monkey1&any_sig), Choice_pref_all_temp(2,monkey1&any_sig),...
            Choice_pref_all_temp(3,monkey1& any_sig),'ko','MarkerSize',14,'MarkerFaceColor',[0.8 0.8 0.8])
        plot3(Choice_pref_all_temp(1,monkey2&any_sig), Choice_pref_all_temp(2,monkey2&any_sig),...
            Choice_pref_all_temp(3,monkey2& any_sig),'k^','MarkerSize',14,'MarkerFaceColor',[0.8 0.8 0.8])
        
        % Orthogonal regression
        [coeff, score, ~, ~, explained,mu] = pca(Choice_pref_all_temp');
        % Direction of the first PC
        dirVect = coeff(:,1); PC1variance = explained(1)/sum(explained);
        
        % Plot the fitted line
        xxx = [min(score(:,1))-.2, max(score(:,1))+.2];
        endpts = [mu + xxx(1)*dirVect' ; mu + xxx(2)*dirVect'];
        plot3(endpts(:,1), endpts(:,2), endpts(:,3), '-', 'color', [0.9290 0.6940 0.1250], 'linew', 3);
        title(['The first PC explains ' num2str(PC1variance) ' variance'])
        
        xlabel('Vestibular'); ylabel('Visual'); zlabel('Combined');
        xticks(-1:0.5:1);         yticks(-1:0.5:1);          zticks(-1:0.5:1); 
        plot3([-1 1], [-1 1],[-1 1], 'k--')
        grid on;
        axis equal; axis square;
        xlim([-1 1]);   ylim([-1 1]);  zlim([-1 1]);
        view(45, 25)
        SetFigure();
        
        
        
        %% ===  Mean raw delta firing (Comb - visual) VS (Comb - vest) %HH20170719 ===
        % %%%%%%%%%% annotated by ZZ 20210111
        %         set(figure(figN),'name','Cdiv(Comb - visual) VS Cdiv(Comb - vest)','pos',[17 514 1151 449]);
        %
        %         % Abs() or not?
        %         mean_raw_delta_firing = mean(PSTH_all_raw_PrefminusNull{1}(:,0<=rate_ts{1} & rate_ts{1}<=1500,:),2);
        %         mean_raw_delta_firing_comb_minus_vest = abs(mean_raw_delta_firing(:,3)) - abs(mean_raw_delta_firing(:,1));
        %         mean_raw_delta_firing_comb_minus_vis = abs(mean_raw_delta_firing(:,3)) - abs(mean_raw_delta_firing(:,2));
        %
        %         h = LinearCorrelation({
        %             (mean_raw_delta_firing_comb_minus_vest(monkey1))
        % %             (mean_raw_delta_firing_comb_minus_vest(monkey2)) ;
        %             },...
        %             {
        %             (mean_raw_delta_firing_comb_minus_vis(monkey1))
        % %             (mean_raw_delta_firing_comb_minus_vis(monkey2)) ;
        %             },...
        %             'CombinedIndex',[],'PlotCombinedOnly', 0, ...
        %             'Xlabel','Combined - Vest (raw PSTH)','Ylabel','Combined - Visual',...
        %             'FaceColors',{'k'},'Markers',{'o','^'},...
        %             'LineStyles',{'k:','k:','k-'},'MarkerSize',marker_size,...
        %             'figN',figN,... 'XHist',20,'YHist',20,'XHistStyle','stacked','YHistStyle','stacked',
        %             'SameScale',1,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        %
        %         delete([h.group(1:2).line]);
        %         plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');         SetFigure(20);
        %         axis([-1 1 -1 1]); axis square;
        %
        %         % Annotate tcells
        %         plot((mean_raw_delta_firing_comb_minus_vest(select_tcells(select_bottom_line))),...
        %             (mean_raw_delta_firing_comb_minus_vis(select_tcells(select_bottom_line))),...
        %             '+','markersize',tcell_cross_size,'color',colors(2,:),'linew',2);
        %         %             plot((Choice_pref_all_temp(2,select_tcells(select_bottom_line),tt)),(Choice_pref_all_temp(3,select_tcells(select_bottom_line),tt)),...
        %         %                 '+','markersize',tcell_cross_size,'color','k','linew',2);
        %
        %         % Show individual cell selected from the figure. HH20150424
        %         h_line = plot([mean_raw_delta_firing_comb_minus_vest],[mean_raw_delta_firing_comb_minus_vis],'visible','off'); hold on;
        %         set([gca h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, [select_cpref_mpref],1});
        %
        %         clear Choice_pref_all_temp;
        %         axis tight;
        
    end

    function f3p2p3(debug)      %  Choice pref (pre vs post)
        if debug
            dbstack;
            keyboard;
        end
        
        %% ---------------- Choice pref (pre vs post) -------------------
        
        set(figure(figN),'name','CDiv (pre-sac) vs. CDiv (post-sac)','pos',[17 514 1151 449]);
        
        k = 3;
        
        cpref_sig_pre = Choice_pref_p_value_all(k,:,1) < 0.01;
        cpref_sig_post = Choice_pref_p_value_all(k,:,3) < 0.01;
        
        h = LinearCorrelation({
            (Choice_pref_all(k,~cpref_sig_pre & ~cpref_sig_post,1)) ;
            (Choice_pref_all(k, xor(cpref_sig_pre , cpref_sig_post),1)) ;
            (Choice_pref_all(k,cpref_sig_pre & cpref_sig_post,1)) ;
            },...
            {
            (Choice_pref_all(k,~cpref_sig_pre & ~cpref_sig_post,2)) ;
            (Choice_pref_all(k, xor(cpref_sig_pre , cpref_sig_post),2)) ;
            (Choice_pref_all(k,cpref_sig_pre & cpref_sig_post,2)) ;
            
            },...
            'CombinedIndex',[7],...
            'Xlabel','Choice preference: decision formation','Ylabel','Choice preference: movement',...
            'FaceColors',{'none',condition_colors(k,:)*0.2 + [0.8 0.8 0.8],condition_colors(k,:)},...
            'EdgeColors',{'k','k','k'},'Markers',{'o'},...
            'LineStyles',{'k:','k:','k:','k-'},'MarkerSize',marker_size,...
            'figN',figN, 'XHist',20,'YHist',20,'XHistStyle','stacked','YHistStyle','stacked',...
            'SameScale',1,'Method','Spearman','FittingMethod',2); figN = figN + 1;         SetFigure(20);
        
        
        delete([h.group(1:3).line]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');
        axis([-1 1 -1 1]); set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1); axis square;
        
        %         % Annotate tcells
        %         plot(Choice_pref_all(k,select_tcells(select_bottom_line),1),Choice_pref_all(k,select_tcells(select_bottom_line),2),...
        %             '+','markersize',tcell_cross_size,'color','k','linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot(Choice_pref_all(k,:,1),Choice_pref_all(k,:,2),'visible','off'); hold on;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_cpref_mpref});
        
        %     for i = 1:3
        %         set(h.group(i).dots,'color',colors(k,:));
        %     end
        
    end

    function f3p3(debug)      % Correlations 3. Psychophysics v.s. CD
        if debug
            dbstack;
            keyboard;
        end
        %%
        % T-cell is important here because Polo's behavior was getting worse while
        % I was getting better at finding T-cells, so...
        select_psycho = select_bottom_line & select_tcells; %& ~[zeros(80,1);ones(138-80,1)];
        
        j = 1;
        
        % Enhancement of cDiv in combined condition
        % enhance_cdiv = max(ChoiceDiv_ModDiffer{1}(:,0 <= rate_ts{j} & rate_ts{j} <= 1500,2),[],2); % Comb - vis
        
        t_begin = 750; t_end = 1050;
%         enhance_cdiv = nanmean(ChoiceDiv_ModDiffer{1}(:,700 <= rate_ts{j} & rate_ts{j} <= 800, 3 ),2); % Comb - vis
        
        h = LinearCorrelation({
            Psy_pred_ratio(select_psycho)
            Psy_pred_ratio(select_psycho)
            Psy_pred_ratio(select_psycho)},...
            {
            nanmean(ChoiceDiv_ModDiffer{1}(select_psycho,t_begin <= rate_ts{j} & rate_ts{j} <= t_end, 1 ),2);
            nanmean(ChoiceDiv_ModDiffer{1}(select_psycho,t_begin <= rate_ts{j} & rate_ts{j} <= t_end, 2 ),2);
            nanmean(ChoiceDiv_ModDiffer{1}(select_psycho,t_begin <= rate_ts{j} & rate_ts{j} <= t_end, 3 ),2);
            },...
            'FaceColors',{condition_colors(1,:),condition_colors(2,:),'k'},'Markers',{'o'},...
            'LineStyles',{'b-','r-','k-'},'MarkerSize',marker_size,...
            'Ylabel',sprintf('\\Delta CDiv (%g - %g ms, 3-1, 3-2, 1-2)',t_begin,t_end),'Xlabel','Psycho prediction ratio',...
            'EdgeColors',{'k','k','k'}, 'MarkerSize',12,...
            'figN',figN,'XHist',15,'YHist',15,'logx',1,...
            'XHistStyle','stacked','YHistStyle','stacked','Method','Pearson','FittingMethod',2); figN = figN + 1;
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');         SetFigure(20);
        
    end

select_group_position; 
    function f3p4(debug)      % Correlations 4. Cell position distribution across longitudinal axes
        if debug
            dbstack;
            keyboard;
        end
        
        position = select_group_position;
        
        % For specific monkey and hemisphere
        to_plot = {monkey_included_for_analysis;1};
        
        if length(monkey_included_for_analysis ) > 1
            select_mon_hems_ind = (any(position(:,1)==monkey_included_for_analysis,2) & position(:,2)==1); % I only select left hemisphere
        elseif length(monkey_included_for_analysis) ==1
            select_mon_hems_ind = (position(:,1)==monkey_included_for_analysis & position(:,2)==1); % I only select left hemisphere
        else
            error('Please choose a monkey !!');
        end
        
        tt = 1;  % Stimulus period
        cpref_sig = (Choice_pref_p_value_all(:,:,tt) < 0.01)';
        mpref_sig = (Modality_pref_p_value_all(:,:,tt) < 0.01)';
        
        Modality_pref_all_stim = (Modality_pref_all(:,:,tt))';     % Modality preference during stimulus period
        Choice_pref_all_stim = (Choice_pref_all(:,:,tt))';
        
        %% For Cells with significant Modality Preference
        
        visual_pref_ind = mpref_sig(:,1) & Modality_pref_all_stim(:,1) >0;
        vestibular_pref_ind = mpref_sig(:,1) & Modality_pref_all_stim(:,1) <0;
        null_pref_ind = ~mpref_sig(:,1);
        
        vispref = position(visual_pref_ind & select_mon_hems_ind,:);
        vestpref = position(vestibular_pref_ind & select_mon_hems_ind,:);
        nullpref = position(null_pref_ind & select_mon_hems_ind,:);
        
        % Mann-Whitney U-test
        AP_sig = ranksum(vispref(:,3),vestpref(:,3));
        depth_sig = ranksum(vispref(:,4),vestpref(:,4));
        
        %== Plotting
        nbins = 10;
        % Manually set bins for AP axis
        min_AP = min(position(:,3));
        max_AP = max(position(:,3));
        xedges_AP = linspace(floor(min_AP),ceil(max_AP), nbins);
        xbins_AP = (xedges_AP(1:end-1)+xedges_AP(2:end))/2;
        
        vispref_num = hist(vispref(:,3),xbins_AP);
        vestpref_num = hist(vestpref(:,3),xbins_AP);
        
        set(figure(figN), 'pos',[50 50 1600 800], 'Name','Distribution of Cell with Modality Preference'); clf;  figN=figN+1;
        hs = subplot(1,2,1); hold on;
        h_bar = bar(xbins_AP,vispref_num, 1,'histc');
        set(h_bar, 'EdgeColor','k',  'FaceColor','r');
        h_bar = bar(xbins_AP,-vestpref_num, 1, 'histc');
        set(h_bar, 'EdgeColor','k',  'FaceColor','b');
        
        set(gca, 'xDir', 'reverse');    % Left -- Anterior
        title('Along AP Axis')
        xlabel({'AC','( mm )'});
        %         if AP_sig <= 0.05
        text(max(xbins_AP),max(vispref_num),num2str(roundn(AP_sig,-4)));
        %         end
        
        
        % Manually set bins for depth (DV axis)
        min_depth = min(position(:,4));
        max_depth = max(position(:,4));
        xedges_depth = (floor(min_depth/500):ceil(max_depth/500))*500;
        xbins_depth = ((xedges_depth(1:end-1))+(xedges_depth(2:end)))/2;
        
        vispref_num = hist(vispref(:,4),xbins_depth);
        vestpref_num = hist(vestpref(:,4),xbins_depth);
        
        hs = subplot(1,2,2); hold on;
        h_bar = bar(xbins_depth,vispref_num, 'histc');
        set(h_bar, 'EdgeColor','k',  'FaceColor','r');
        h_bar = bar(xbins_depth,-vestpref_num, 'histc');
        set(h_bar, 'EdgeColor','k',  'FaceColor','b');
        
        legend({'Visual'; 'Vestibular'}, 'Position',[0.9 0.9 0.1 0.1],'AutoUpdate','off');
        title('Along DV Axis');
        xlabel({'Depth into CD', '( \mum)'});
        %         if depth_sig <= 0.05
        text(0,max(vispref_num),num2str(roundn(depth_sig,-4)));
        %         end
        
        set(gca, 'view', [-90,-90]);
        
        SetFigure();
        
        %% For cells with significance Choice Divergence
        
        any_cell = any(cpref_sig')';
        anycell_position = position(any_cell, :);
        noncell_position = position(~any_cell, :);
        
        anycell_AP_counts = hist(anycell_position(:,3),xbins_AP);
        anycell_depth_counts = hist(anycell_position(:,4),xbins_depth);
        noncell_AP_counts = hist(noncell_position(:,3),xbins_AP);
        noncell_depth_counts = hist(noncell_position(:,4),xbins_depth);
        
        % Mann-Whitney U-test
        AP_sig = ranksum(anycell_position(:,3),noncell_position(:,3));
        depth_sig = ranksum(anycell_position(:,4),noncell_position(:,4));
        
        set(figure(figN), 'pos',[70 50 1600 800], 'Name','Distribution of Cell with Choice Preference'); clf;  figN=figN+1;
        hs = subplot(1,2,1); hold on;
        h_bar = bar(xbins_AP,anycell_AP_counts, 'histc');
        set(h_bar, 'EdgeColor','k',  'FaceColor','k');
        h_bar = bar(xbins_AP,-noncell_AP_counts, 'histc');
        set(h_bar, 'EdgeColor','k',  'FaceColor',[0.5,0.5,0.5]);
        
        set(gca, 'xDir', 'reverse');    % Left -- Anterior
        title('Along AP Axis')
        xlabel({'AC','( mm )'});
        %         if AP_sig <= 0.05
        text(max(xbins_AP),max(anycell_AP_counts),num2str(roundn(AP_sig,-4)));
        %         end
        
        hs = subplot(1,2,2); hold on;
        h_bar = bar(xbins_depth,anycell_depth_counts, 'histc');
        set(h_bar, 'EdgeColor','k',  'FaceColor','k');
        h_bar = bar(xbins_depth,-noncell_depth_counts, 'histc');
        set(h_bar, 'EdgeColor','k',  'FaceColor',[0.5,0.5,0.5]);
        
        legend({'Any-Cell'; 'Non-Cell'},'Position',[0.9 0.9 0.1 0.1],'AutoUpdate','off');
        title('Along DV Axis');
        xlabel({'Depth into CD', '( \mum)'});
        %         if depth_sig <= 0.05
        text(0,max(anycell_depth_counts),num2str(roundn(depth_sig,-4)));
        %         end
        
        set(gca, 'view', [-90,-90]);
        
        SetFigure();
        
        
    end

ChoiceDiv_All_perm_select; 
    function f4p0(debug)      % Pack PCA_A
        if debug
            dbstack;
            keyboard;
        end
        %%
        % Pack all data into matrix A: [memsacDDI, vest_div, vis_div, comb_div]
        
        % Now it becomes A: [memsacDDI, vest_div, vis_div, comb_div, 2-1 mod_div, 3-1 mod_div, 3-2 mod_div]
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        sort_time_interval1 = [find(rate_ts{j_PCA_A} > time_markers{j_PCA_A}(1,1), 1)  find(rate_ts{j_PCA_A} >= time_markers{j_PCA_A}(1,2),1)-1]; % Stim on to saccade on % changed to stim on to stim off
        sort_time_interval2 = [ find(rate_ts{2}>=min(rate_ts{2}),1) find(rate_ts{2} > time_markers{2}(1,1),1)-1]; % Post sac % changed to pre sac
        sort_time_interval3 = [1  length(CP_ts{j_PCA_A})];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Exclude any cells that have NaNs
        % for k = 1:3
        %     selectCells = selectCells & ~isnan(ChoiceDiv_All(:,1,k));
        % end
        
        
        % --- Sort A according to different methods for visualizing ---
        
        A_CP = reshape(CP{j_PCA_A}(select_bottom_line,:,:),sum(select_bottom_line),[])-0.5;
        
        A_choicediv= reshape(ChoiceDiv_All{j_PCA_A}(select_bottom_line,:,:),sum(select_bottom_line),[]);
        
        % A_memSac = xls_num{1}(selectCells,header.DDI_LS:header.DDI_post) - 0.4; % All memSac data
        % LS, M, Pre, Co, Post
        A_memSac = group_MemSac_DDI(select_bottom_line,2:end);  % Read memsac from .mat file. HH20150413
        A_memSac = (A_memSac-0.5)*2; % Map to -1 ~ 1 for better visualization
        
        % Now I add mod_div. HH20150415
        % A_moddiv = reshape(ModDiv_All{j_PCA_A}(select_bottom_line,:,:),sum(select_bottom_line),[]);
        A_moddiv = reshape(ModDiv_All{j_PCA_A}(select_bottom_line,:,1),sum(select_bottom_line),[]); % Only vis-vest
        
        PCA_A = [A_memSac A_choicediv A_moddiv]; % A_CP];
        
    end
    function f4p1(debug)      % 1. Hot-gram
        if debug
            dbstack;
            keyboard;
        end
        
        % If data not packed
        if isempty(PCA_A)
            f4p0(0);
        end
        
        sort_method = {
            [2 2];
            %             [1 4];
            %     [3 4];
            
            5 + 2*length(rate_ts{j_PCA_A}) + sort_time_interval1;  % sorted by mean CD during stim period in combined condition
            
            5 + 3* length(rate_ts{j_PCA_A}) + sort_time_interval1; %
            };
        
        for sort_ind = length(sort_method):-1:1
            
            sort_begin = sort_method{sort_ind}(1);
            sort_end = sort_method{sort_ind}(end);
            
            mean_for_sort = mean(PCA_A(:,sort_begin:sort_end),2);
            mean_for_sort(isnan(mean_for_sort)) = -inf;  % Nan goes last
            
            [~,sort_order] = sort(mean_for_sort,'descend');
            
            % Enlarge mem-sac part for clarity
            A_forplot = [reshape(repmat(A_memSac,enlarge_factor,1),size(A_memSac,1),[]) A_choicediv A_moddiv A_moddiv(:,end)];% A_CP A_CP(:,end)]; % To ensure the last value to be shown
            A_forplot = A_forplot(sort_order,:);
            
            A_forplot = [A_forplot; A_forplot(end,:)]; % To ensure the last cell to be shown
            
            A_forplot(isnan(A_forplot)) = -0.5; % This workaround is to avoid that cells which are next the NaNs cannot be displayed.
            
            % [nn, tt] = meshgrid(1:size(A_forplot,1), 1:size(A_forplot,2));
            
            set(figure(2099+figN),'name',['CDiv and MDiv, hot-gram, j_PCA_A = ' num2str(j_PCA_A)],...
                'unit','norm','pos',[0.581547619047619 0.177142857142857 0.414880952380952 0.738095238095238]);
            clf; figN = figN+1;
            
            % h1 = surf(tt,nn,A_forplot','Edgecolor','none');
            h1 = imagesc(A_forplot); colormap jet
            % axis xy;
            
            % Annotate time separations
            hold on;
            CD_begin_at = enlarge_factor * 5 + 1;
            
            % Stim on / stim off / sac on
            for ttt = 1:2
                for tt = 0: 2% 5
                    plot(CD_begin_at + tt*length(rate_ts{j_PCA_A}) + find(rate_ts{j_PCA_A} >= time_markers{j_PCA_A}(1,ttt),1) * [1 1],...
                        [1 size(A_forplot,1)],'k','linesty',marker_for_time_markers{j_PCA_A}{ttt},'linewid',1.5);
                end
            end
            
            plot3([CD_begin_at CD_begin_at],[1 size(A_forplot,1)],[1 1],'k','linewid',3);
            for tt = 1: 4 % 6
                plot(CD_begin_at + tt*[length(rate_ts{j_PCA_A})  length(rate_ts{j_PCA_A})],[1 size(A_forplot,1)],'k','linewid',3);
            end
            
            % Annotate tcells
            cell_loc = select_tcells(select_bottom_line);
            h2 = plot(0,find(cell_loc(sort_order))+.5,'+','markersize',5,'color','k','linew',1.5);
            
            % Annotate sort methods
            temp_beg = @(x)((x<=5)*((x-1)*enlarge_factor) + 1 + (x>5)*(CD_begin_at+(x-5)-1));
            temp_end = @(x)((x<=5)*((x)*enlarge_factor) + 1 + (x>5)*(CD_begin_at+(x-5)));
            
            sort_begin_forplot = temp_beg(sort_begin);
            sort_end_forplot = temp_end(sort_end);
            
            plot([sort_begin_forplot sort_end_forplot],[size(A_forplot,1)+1 size(A_forplot,1)+1],'color',condition_colors(2,:),'linewid',5);
            xlabel('Temporal features');
            ylabel('Cell Number');
            
            set(gca,{'xtick','ytick','ztick'},{[],10:10:size(A_forplot,1),[]},'linewidth',0.00000001);
            axis tight; colorbar; ylim([1 size(PCA_A,1)+2.5]); SetFigure();
            
            % Show individual cell selected from the figure. HH20150424
            [~,where] = sort(sort_order);
            h_line = plot(zeros(1,length(where)),where+.5,'visible','off'); hold on;
            set([gca; h1; h2(:); h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_for_PCA_B});
            
            
            % Only plot significant CD distribution to visualize heterogenity
            % Added by ZZ 20220309
            if sort_ind ==2
                set(figure(220309),'name',['CDiv hot-gram, j_PCA_A = ' num2str(j_PCA_A)],...
                    'unit','norm','pos',[0.11547619047619 0.177142857142857 0.714880952380952 0.738095238095238]);clf;
                for c = 1: stim_type
                    subplot(1,3,c); hold on;
                    
                    sig_ind = ChoiceDiv_All_perm_select{j_PCA_A}.p <= 0.01;
                    
                    h1 = imagesc(sig_ind(sort_order,:,c));
                    colormap parula;
                    
                    for ttt = 1:2
                        plot(find(rate_ts{j_PCA_A} >= time_markers{j_PCA_A}(1,ttt),1) * [1 1],...
                            [1 size(sig_ind,1)],'k','linesty',marker_for_time_markers{j_PCA_A}{ttt},'linewid',1.5);
                    end
                    
                    % Annotate tcells in all conditions
                    cell_loc = group_ChoicePreference_pvalue(c,select_bottom_line, j_PCA_A) <= 0.01;
                    h2 = plot(size(sig_ind,2)+5, find(cell_loc(sort_order)), 's',...
                        'MarkerSize', 4, 'Color', condition_colors(c,:), 'MarkerFaceColor',condition_colors(c,:));
                    
                    xlabel('Time to stimulus onset (ms)');
                    ylabel('Cell Number');
                    axis tight;
                    ylim([1 size(sig_ind,1)]); 
                    xticks([]);
                    
                    % Gaussian vel
                    % Only extract 150 time point
                    extract_Gauss_vel = Gauss_vel(ceil(linspace(1,length(Gauss_vel),150)),2);
                    extract_time_bin = (1:150) + find(rate_ts{j_PCA_A} >= time_markers{j_PCA_A}(1,1),1);
                    
                    plot(extract_time_bin',extract_Gauss_vel*max(ylim)/5,'--','linew',3,'color',[0.6 0.6 0.6]);
                    %         axis tight
                    
                end
                
                SetFigure();
            end
            
        end
        
    end
    function f4p2(debug)      % 2. PCA_A (Eigenfeature)
        if debug;  dbstack;  keyboard;  end
        
        % If data not packed
        if isempty(PCA_A)
            f4p0(debug);
        end
        
        % Use the original A instead of A_forplot
        dataPCA = PCA_A;
        
        dataPCA(sum(isnan(dataPCA),2)>0,:) = [];
        
        % [~,sort_order] = sort(mean(dataPCA(:,3:4),2)); % Sort according to mem-sac
        [~,sort_order] = sort(mean(dataPCA(:, 3:4 + 2*length(rate_ts{j_PCA_A}) + sort_time_interval1),2));  % sorted by CD in combined condition
        dataPCA = dataPCA(sort_order,:);
        
        % Do PCA
        [PCA_A_PC, ~] = pca(dataPCA);
        
        % The first three eigenvectors that have the first three largest
        % eigenvalues.
        PC1 = PCA_A_PC(:,1);
        PC2 = PCA_A_PC(:,2);
        % PC3 = sortedEigVectors(:,3);
        
        % Projecting the raw data onto the first three eigenvectors
        projPC1 = dataPCA * PC1;
        projPC2 = dataPCA * PC2;
        % projPC3 = dataPCA(:,2:end) * PC3;
        
        set(figure(2099+figN),'name',['Scatter of Eigenfeature, j_PCA_A = ' num2str(j_PCA_A)],'pos',[41 445 543 451]); clf; figN = figN+1;
        
        colorOrder = jet(length(projPC1));
        
        scatter(projPC1, projPC2, 150, colorOrder,'fill');
        grid on; hold on;
        xlabel('PC1'); ylabel('PC2');
        SetFigure();
        
        set(figure(2099+figN),'pos',[602 446 978 450],'name',['Eigenfeature, j_PCA_A = ' num2str(j_PCA_A)]); clf; figN = figN+1;
        
        for ss = 1:2
            subplot(2,1,ss);
            
            eval(['plot(PC' num2str(ss) ',''lineW'',2)']);
            
            hold on;
            
            plot([5 5],ylim,'k','linewid',1);
            plot(xlim,[0 0],'k--');
            
            % Stim on / stim off / sac on
            for ttt = 1:2
                for tt = 0:3
                    plot(5 + tt*length(rate_ts{j_PCA_A}) + find(rate_ts{j_PCA_A} >= time_markers{j_PCA_A}(1,ttt),1) * [1 1],...
                        ylim,'k','linesty',marker_for_time_markers{j_PCA_A}{ttt},'linewid',1.5);
                end
            end
            
            % End
            plot([5 5],ylim,'k','linewid',3);
            for tt = 1:4
                plot(5 + tt*[length(rate_ts{j_PCA_A})  length(rate_ts{j_PCA_A})],ylim,'k','linewid',3);
            end
            
            axis tight; axis off
        end
        
        SetFigure();
        
    end

    function f5p0(debug)      % Do PCA_B analysis
        if debug;  dbstack;  keyboard;  end
        
        
        %% Do PCA_B analysis
        % Pack all data (PSTH) into matrix B: [vest_pref, vest_null, vis_pref, vis_null, comb_pref, comb_null]
        % Only decicion period is included
        
        find_PCA_B = find(select_for_PCA_B);
        
        PCA_B = nan(sum(select_for_PCA_B),6 * sum(PCA_B_time_range));
        
        for i = 1:sum(select_for_PCA_B)
            raw_PSTH = group_result(find_PCA_B(i)).mat_raw_PSTH.PSTH{j_PCA_B,ALL_CorrectCHOICE,1}.ys;
            if size(raw_PSTH,1)==6
                PCA_B(i,:) = reshape(raw_PSTH(:,PCA_B_time_range)',[],1)';
            end
        end
        
        PCA_B(sum(isnan(PCA_B),2) > 0,:) = [];
        
        % Do PCA
        [weights_PCA_B_PC, score, latent, ~, PCA_B_explained] = pca(PCA_B');
        
        % Projecting the raw data onto the first several eigenvectors
        for dim = 1:denoised_dim
            PCA_B_projPC{dim} = (reshape(weights_PCA_B_PC(:,dim)' * PCA_B,[],6))';
            %     projPC{dim} = reshape(score(:,dim),[],6)';
        end
        
    end
    function f5p1(debug)      % 1. Weights and correlations  % HH20150413
        if debug;  dbstack;  keyboard;  end
        
        % If PCA_B has never been done, do it first
        if isempty(PCA_B_projPC)
            f5p0(0);
        end
        
        % ------- 3-d PCA weights ---------
        set(figure(2099+figN),'pos',[13 53 836 702],'name','PCA weights in 3d'); clf; figN = figN+1; hold on;
        
        tcell_in_bottom_line = select_tcells(select_for_PCA_B);
        plot3(weights_PCA_B_PC(tcell_in_bottom_line,1),weights_PCA_B_PC(tcell_in_bottom_line,2),weights_PCA_B_PC(tcell_in_bottom_line,3),'r+','markersize',15,'linew',1.5);
        
        h_line = plot3(weights_PCA_B_PC(:,1),weights_PCA_B_PC(:,2),weights_PCA_B_PC(:,3),'ko','markersize',13,'linew',1.5);
        grid off;
        SetFigure(); xlabel('Contribution to Eigen-neuron 1 (choice)'); ylabel('Contribution to Eigen-neuron 2 (modality)'); zlabel('PC3'); grid on; axis tight;
        
        % Show individual cell selected from the figure. HH20150424
        set([gca h_line],'ButtonDownFcn',{@Show_individual_cell,h_line,select_for_PCA_B});
        
        % -------- Weights vs cell properties ----------
        Weights_Property_Correlation(weights_PCA_B_PC(:,[1 3]),...
            {'Contribution to Eigen-neuron 1','Contribution to Eigen-neuron 3'},select_for_PCA_B);
        
    end
    function f5p2(debug)      % 2. 1-D Trajectory
        if debug;  dbstack;  keyboard;  end
        
        % If PCA_B has neven been done, do it first
        if isempty(PCA_B_projPC)
            f5p0(0);
        end
        
        % ============  1. Variance explained =================
        set(figure(2099+figN),'pos',[936 491 733 463],'name',['Variance explained, j_PCA_B = ' num2str(j_PCA_B)]); clf; figN = figN+1; hold on;
        plot((1:length(PCA_B_explained))', cumsum(PCA_B_explained),'o-','markersize',8,'linew',1.5);
        plot((1:denoised_dim)',cumsum(PCA_B_explained(1:denoised_dim)),'ro-','markersize',8,'linew',1.5,'markerfacecol',condition_colors(2,:));
        plot([0 1],[0 PCA_B_explained(1)],'r-','linew',1.5);
        plot(xlim,[1 1]*sum(PCA_B_explained(1:denoised_dim)),'r--');
        text(denoised_dim,sum(PCA_B_explained(1:denoised_dim))*0.9,[num2str(sum(PCA_B_explained(1:denoised_dim))) '%'],'color',condition_colors(2,:));
        SetFigure(); xlabel('Num of principal components'); ylabel('Explained variability (%)'); ylim([0 100]);
        
        % ============  2. 1-D Trajectory of Eigen-neurons =================
        
        set(figure(2099+figN),'pos',[462 67 1207 889],'name',['Population Dynamics, j_PCA_B = ' num2str(j_PCA_B)]); clf; figN = figN+1; hold on;
        ds = [1:8];
        
        for d = 1:length(ds)
            ranges(d) = range(PCA_B_projPC{ds(d)}(:));
        end
        
        % Plotting
        for d = 1:length(ds)
            % Normalize to [-1,1]
            gain = max(ranges) / (2 * 0.8) ;
            offset = mean(PCA_B_projPC{ds(d)}(:));
            norm_proj_PC_this = (PCA_B_projPC{ds(d)}-offset)/gain;
            
            h = subplot(fix(sqrt(length(ds))),ceil(length(ds)/fix(sqrt(length(ds)))),d);
            SeriesComparison(shiftdim(norm_proj_PC_this',-1),PCA_B_times,'Colors',{condition_colors(1,:),condition_colors(1,:),condition_colors(2,:),condition_colors(2,:),[0 0.8 0.4],[0 0.8 0.4]},'LineStyles',{'-','--'},'axes',h);
            axis tight; ylim([-1 1]); xlim([-200 1860])
            for tt = 1:2
                plot([1 1] * time_markers{j_PCA_B}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j_PCA_B}{tt},'linew',1.5);
            end
            set(gca,'ytick',[-1 0 1]);
            xlabel([]); ylabel('Amplitude (a.u.)'); title(['Eigen-neuron ' num2str(ds(d))]); legend off;
        end
        % set(get(gcf,'children'),'ylim',[min_ylim max_ylim]);
        SetFigure(15);
        
    end
    function f5p3(debug)      % 3. 2-D Trajectory  % HH20150413
        if debug;  dbstack;  keyboard;  end
        
        % If PCA_B has neven been done, do it first
        if isempty(PCA_B_projPC)
            f5p0(0);
        end
        
        % ======== 2D ========= %  % Increase to 3d
        set(figure(2099+figN),'pos',[18 170 898 786],'name',['Population Dynamics, j_PCA_B = ' num2str(j_PCA_B)]); clf; figN = figN+1; hold on;
        
        which_two_dimension = [1,2,3];
        %         which_two_dimension = [2,3];
        
        % -- For plotting time markers
        plotInt = 150; % ms
        plotPerTimeBin = fix(plotInt/(PCA_B_times(2)-PCA_B_times(1))); % Should be 100/10 = 10
        plotMinInd = fix(-min(PCA_B_times)/plotInt)*plotPerTimeBin;
        plotInd = plotMinInd : plotPerTimeBin : length(PCA_B_times);
        %%
        cla
        for k = 1:3
            
            % Time markers
            start_time = 1; % Start point
            
            h_pref( k) = plot3(PCA_B_projPC{which_two_dimension(1)}((k-1)*2+1,start_time),...
                PCA_B_projPC{which_two_dimension(2)}((k-1)*2+1,start_time),...
                PCA_B_projPC{which_two_dimension(3)}((k-1)*2+1,start_time),...
                'o','color',condition_colors(k,:),'markersize',20,'markerfacecol',condition_colors(k,:));
            % Null
            h_null(k) = plot3(PCA_B_projPC{which_two_dimension(1)}(k*2,start_time),...
                PCA_B_projPC{which_two_dimension(2)}(k*2,start_time),...
                PCA_B_projPC{which_two_dimension(3)}(k*2,start_time),...
                'o','color',condition_colors(k,:),'markersize',20,'linew',3);
            
            
            plot3(PCA_B_projPC{which_two_dimension(1)}((k-1)*2+1,:),...
                PCA_B_projPC{which_two_dimension(2)}((k-1)*2+1,:),...
                PCA_B_projPC{which_two_dimension(3)}((k-1)*2+1,:),...
                '-','color',condition_colors(k,:),'linew',3);
            % Null
            plot3(PCA_B_projPC{which_two_dimension(1)}(k*2,:),...
                PCA_B_projPC{which_two_dimension(2)}(k*2,:),...
                PCA_B_projPC{which_two_dimension(3)}(k*2,:),...
                '--','color',condition_colors(k,:),'linew',3);
            
            
            % -- Time markers
            colorsHsv = repmat(rgb2hsv(condition_colors(k,:)),length(plotInd),1);
            colorsHsv(:,2) = linspace(0.2,1,length(plotInd));
            
            if k == 3
                colorsHsv(:,3) = linspace(.9,colorsHsv(1,3),length(plotInd));
            end
            
            colorsRGB = hsv2rgb(colorsHsv);
            
            for pp = 1:length(plotInd)
                plot3(PCA_B_projPC{which_two_dimension(1)}((k-1)*2+1,plotInd(pp)),...
                    PCA_B_projPC{which_two_dimension(2)}((k-1)*2+1,plotInd(pp)),...
                    PCA_B_projPC{which_two_dimension(3)}((k-1)*2+1,plotInd(pp)),...
                    'o','color',colorsRGB(pp,:),'markerfacecol',colorsRGB(pp,:),'linew',0.1,'markersize',15);
                
                plot3(PCA_B_projPC{which_two_dimension(1)}(k*2,plotInd(pp)),...
                    PCA_B_projPC{which_two_dimension(2)}(k*2,plotInd(pp)),...
                    PCA_B_projPC{which_two_dimension(3)}(k*2,plotInd(pp)),...
                    'o','color',colorsRGB(pp,:),'markerfacecol','none','linew',2,'markersize',15);
                
            end
            
        end
        %%
        %         axis tight;
        grid off;
        axis off
        % xlabel('PC1'); ylabel('PC2');
        
        % xlims = xlim; ylims = ylim;
        % h_text = text(xlims(1),ylims(2),'');
        
        % ======= Euclidean distance =========
        h_timeaxis = axes('pos', [0.026 0.616 0.25 0.3] ,'color','none');  hold on
        
        prefs = [1 3 5];
        nulls = [2 4 6];
        
        % distance = sqrt((PCA_B_projPC{1}(prefs,:) - PCA_B_projPC{1}(nulls,:)).^2 + (PCA_B_projPC{2}(prefs,:) - PCA_B_projPC{2}(nulls,:)).^2 + (PCA_B_projPC{3}(prefs,:) - PCA_B_projPC{3}(nulls,:)).^2);
        % distance = distance';
        
        distance = sum(reshape(cell2mat((cellfun(@(x)(x(prefs,:)-x(nulls,:)).^2',PCA_B_projPC,'uni',0))),[],3,denoised_dim),3);
        
        % Normalize distance
        % distance = (distance - repmat(min(distance,[],1),size(distance,1),1))./repmat(range(distance),size(distance,1),1);
        % set(figure(2099+figN),'pos',[18 355 766 601],'name',['Euclidean distance, j_PCA_B = ' num2str(j_PCA_B)]); clf; figN = figN+1; hold on;
        
        for k = 1:3
            plot(repmat(PCA_B_times',1,3),distance(:,k),'color',condition_colors(k,:),'linew',3);
        end
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{j_PCA_B}(1),Gauss_vel(:,2)*max(ylim)/4,'--','linew',3,'color',[0.6 0.6 0.6]);
        axis tight
        
        %         plot(repmat(PCA_B_times',1,3),sum(distance(:,[1 2]),2),'k--','linew',2);
        %         ylim([0 max(distance(:))]);
        
        for tt = 1:2
            plot([1 1] * time_markers{j_PCA_B}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j_PCA_B}{tt},'linew',2);
        end
        
        text(mean(xlim)/5,max(ylim)*1.05, sprintf('Euclidean distance\nin %g-d space \n(N = %g)', denoised_dim, sum(select_for_PCA_B)));
        
        SetFigure(15);
        
        % =========== Animation ===========
        
        h_timeline = plot([PCA_B_times(1) PCA_B_times(1)],ylim,'m','linew',3);
        
        set(h_timeaxis,'xlim',[PCA_B_times(1)  PCA_B_times(end)],'ytick',[],'ycolor','w');
        
        uicontrol('Style','pushbutton','Unit','norm','Position',[0.032 0.022 0.085 0.054],...
            'String','Play','FontSize',13,...
            'Callback',{@Population_dynamics,0}); % Play
        uicontrol('Style','pushbutton','Unit','norm','Position',[0.132 0.022 0.085 0.054],...
            'String','Record','FontSize',13,...
            'Callback',{@Population_dynamics,10}); % Record
        uicontrol('Style','slider','Unit','norm','Position',[0.026 0.53 0.2 0.03],...
            'Min',1,'Max',length(PCA_B_times),'Value',1,...
            'Callback',{@Population_dynamics,999});
        
        % ======= Momentum ========
        figure(532);  clf;
        
        plot(repmat(PCA_B_times(1:end-1)',1,3), diff(distance),'linew',3); hold on;
        axis([-50 1800 -300 500]);
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{j_PCA_B}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',2,'color',[0.6 0.6 0.6]);
        for tt = 1:2
            plot([1 1] * time_markers{j_PCA_B}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j_PCA_B}{tt},'linew',2);
        end
        SetFigure(15);
        
        % Note here I use a double nested function. @HH20150425
        % Play population dynamics. HH20150424
        function Population_dynamics(~,~,flag)
            switch flag
                case {0,10}   % Play & record
                    % Preallocate movie structure.
                    mov(1:size(PCA_B_projPC{1},2)) = struct('cdata', [],...
                        'colormap', []);
                    
                    img = VideoWriter('Caudate_dynamics_j_1.avi', 'Uncompressed AVI');
                    open(img);
                    
                    for t_ind = 1 : size(PCA_B_projPC{1},2)
                        for kk = 1:3
                            % Update
                            set(h_pref(kk),'xdata',PCA_B_projPC{which_two_dimension(1)}((kk-1)*2+1,t_ind),...
                                'ydata',PCA_B_projPC{which_two_dimension(2)}((kk-1)*2+1,t_ind),...
                                'zdata',PCA_B_projPC{which_two_dimension(3)}((kk-1)*2+1,t_ind));
                            set(h_null(kk),'xdata',PCA_B_projPC{which_two_dimension(1)}(kk*2,t_ind),...
                                'ydata',PCA_B_projPC{which_two_dimension(2)}(kk*2,t_ind),...
                                'zdata',PCA_B_projPC{which_two_dimension(3)}(kk*2,t_ind));
                        end
                        
                        %     set(h_text,'string',num2str(fix(rate_ts{j_PCA_B}(tt)/10)*10));
                        set(h_timeline,'xdata',PCA_B_times(t_ind)*[1 1]);
                        
                        if flag
                            mov(t_ind) = getframe(gcf);
                            writeVideo(img, mov(t_ind));
                        else
                            pause(0.01);
                        end
                        
                        drawnow;
                    end
                                        
                case 999  % From slider bar
                    t_ind = round(get(gcbo,'value'));
                    
                    for kk = 1:3
                        % Update
                        set(h_pref(kk),'xdata',PCA_B_projPC{which_two_dimension(1)}((kk-1)*2+1,t_ind),...
                            'ydata',PCA_B_projPC{which_two_dimension(2)}((kk-1)*2+1,t_ind),...
                            'zdata',PCA_B_projPC{which_two_dimension(3)}((kk-1)*2+1,t_ind));
                        set(h_null(kk),'xdata',PCA_B_projPC{which_two_dimension(1)}(kk*2,t_ind),...
                            'ydata',PCA_B_projPC{which_two_dimension(2)}(kk*2,t_ind),...
                            'zdata',PCA_B_projPC{which_two_dimension(3)}(kk*2,t_ind));
                    end
                    
                    %     set(h_text,'string',num2str(fix(rate_ts{j_PCA_B}(tt)/10)*10));
                    set(h_timeline,'xdata',PCA_B_times(t_ind)*[1 1]);
                    
            end
        end
        
    end

    function f51p1(debug)     % 1. PCA over heading and choice for each modality (Dora's method) %HH20160701
        if debug;  dbstack;  keyboard;  end
        
        find_PCA_B = find(select_for_PCA_B);
        
        for k = 1:3
            PCA_B_heading_eachk{k} = PSTH_correct_angles_raw{j_PCA_B}(select_for_PCA_B,PCA_B_time_range,:,k);
            PCA_B_heading_eachk{k} = reshape(PCA_B_heading_eachk{k},size(PCA_B_heading_eachk{k},1),[]); % Cell * (0-pref,0-null,1-pref,1-null, ..., 8-pref, 8-null)
        end
        
        PCA_B_heading = [PCA_B_heading_eachk{1} PCA_B_heading_eachk{2} PCA_B_heading_eachk{3}];
        PCA_B_heading(sum(isnan(PCA_B_heading),2) > 0,:) = [];
        
        % Do PCA
        [weights_PCA_B_heading_PC, score, latent, ~, PCA_B_explained_heading] = pca(PCA_B_heading');
        
        % Manually rotate PC1 and PC2 (they're also orthogonal, so there's no
        % change in PCA results). @HH20150424
        %         theta = 55 / 180 * pi; % Rotate counterclockwise the two PC bases
        %         THETA = [cos(theta) -sin(theta) ; sin(theta) cos(theta)];
        %         weights_PCA_B_PC(:,1:2) = weights_PCA_B_PC(:,1:2) * THETA;
        %         weights_PCA_B_PC(:,2) = -weights_PCA_B_PC(:,2);
        
        % Projecting the raw data onto the first several eigenvectors
        for dim = 1:denoised_dim
            PCA_B_heading_projPC{dim} = (reshape(weights_PCA_B_heading_PC(:,dim)' * PCA_B_heading,[],3*size(PSTH_correct_angles_raw{j_PCA_B},3)))';
            %     projPC{dim} = reshape(score(:,dim),[],6)'
        end
        
        % ======== 2D ========= %
        set(figure(2099+figN),'pos',[18 170 898 786],'name',['Population Dynamics, j_PCA_B = ' num2str(j_PCA_B)]); clf; figN = figN+1; hold on;
        
        which_two_dimension = [2,3,1];  % changed by ZZ 20210511
        %         which_two_dimension = [2,3];
        
        for k = 1:3
            %             k = 1;
            colors_angles{k} = colormap(gray);
            colors_angles{k} = ones(size(colors_angles{k},1),3) - colors_angles{k} .* repmat([1 1 1]-condition_colors(k,:),size(colors_angles{k},1),1);
            colors_angles{k} = colors_angles{k}(round(linspace(20,length(colors_angles{k}),5)),:);
            
            for hh = 1:5 % 5 headings
                %             % Time markers
                start_time = 1; % Start point
                % Pref
                plot3(PCA_B_heading_projPC{which_two_dimension(1)}(2*hh-1+(k-1)*10,start_time),...
                    PCA_B_heading_projPC{which_two_dimension(2)}(2*hh-1+(k-1)*10,start_time),...
                    PCA_B_heading_projPC{which_two_dimension(3)}(2*hh-1+(k-1)*10,start_time),...
                    'o','color',colors_angles{k}(hh,:),'markersize',20,'markerfacecol',colors_angles{k}(hh,:));
                % Null
                plot3(PCA_B_heading_projPC{which_two_dimension(1)}(2*hh+(k-1)*10,start_time),...
                    PCA_B_heading_projPC{which_two_dimension(2)}(2*hh+(k-1)*10,start_time),...
                    PCA_B_heading_projPC{which_two_dimension(3)}(2*hh+(k-1)*10,start_time),...
                    'o','color',colors_angles{k}(hh,:),'markersize',20,'linew',3);
                
                % Pref
                plot3(PCA_B_heading_projPC{which_two_dimension(1)}(2*hh-1+(k-1)*10,:),...
                    PCA_B_heading_projPC{which_two_dimension(2)}(2*hh-1+(k-1)*10,:),...
                    PCA_B_heading_projPC{which_two_dimension(3)}(2*hh-1+(k-1)*10,:),...
                    '-','color',colors_angles{k}(hh,:),'linew',3);
                % Null
                plot3(PCA_B_heading_projPC{which_two_dimension(1)}(2*hh+(k-1)*10,:),...
                    PCA_B_heading_projPC{which_two_dimension(2)}(2*hh+(k-1)*10,:),...
                    PCA_B_heading_projPC{which_two_dimension(3)}(2*hh+(k-1)*10,:),...
                    '--','color',colors_angles{k}(hh,:),'linew',3);
            end
            
            axis tight;  grid off;
            axis off
        end
        
        
        % ============  1. Variance explained =================
        set(figure(2099+figN),'pos',[936 491 733 463],'name',['Variance explained, j_PCA_B = ' num2str(j_PCA_B)]); clf; figN = figN+1; hold on;
        plot((1:length(PCA_B_explained_heading))', cumsum(PCA_B_explained_heading),'o-','markersize',8,'linew',1.5);
        plot((1:denoised_dim)',cumsum(PCA_B_explained_heading(1:denoised_dim)),'ro-','markersize',8,'linew',1.5,'markerfacecol',condition_colors(2,:));
        plot([0 1],[0 PCA_B_explained_heading(1)],'r-','linew',1.5);
        plot(xlim,[1 1]*sum(PCA_B_explained_heading(1:denoised_dim)),'r--');
        text(denoised_dim,sum(PCA_B_explained_heading(1:denoised_dim))*0.9,[num2str(sum(PCA_B_explained_heading(1:denoised_dim))) '%'],'color',condition_colors(2,:));
        SetFigure(); xlabel('Num of principal components'); ylabel('Explained variability (%)'); ylim([0 100]);
        
        % ============  2. 1-D Trajectory of Eigen-neurons =================
        
        set(figure(2099+figN),'pos',[462 67 1207 889],'name',['Population Dynamics, j_PCA_B = ' num2str(j_PCA_B)]); clf; figN = figN+1; hold on;
        ds = [1:8];
        
        for d = 1:length(ds)
            ranges(d) = range(PCA_B_heading_projPC{ds(d)}(:));
        end
        
        % Plotting
        for d = 1:length(ds)
            % Normalize to [-1,1]
            gain = max(ranges) / (2 * 0.8) ;
            offset = mean(PCA_B_heading_projPC{ds(d)}(:));
            norm_proj_PC_this = (PCA_B_heading_projPC{ds(d)}-offset)/gain;
            
            h = subplot(fix(sqrt(length(ds))),ceil(length(ds)/fix(sqrt(length(ds)))),d);
            
            colors_angles_all = [colors_angles{1};colors_angles{2};colors_angles{3}];
            colors_angles_all = reshape(repmat(colors_angles_all,1,2)',3,[])';
            colors_angles_all = mat2cell(colors_angles_all,ones(30,1));
            
            SeriesComparison(shiftdim(norm_proj_PC_this',-1),PCA_B_times,'Colors',...
                colors_angles_all,'LineStyles',{'-','--'},'axes',h);
            axis tight; ylim([-1 1]); xlim([-200 1700])
            for tt = 1:2
                plot([1 1] * time_markers{j_PCA_B}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j_PCA_B}{tt},'linew',1.5);
            end
            set(gca,'ytick',[-1 0 1]);
            xlabel([]); ylabel('Amplitude (a.u.)'); title(['Eigen-neuron ' num2str(ds(d))]); legend off;
        end
        % set(get(gcf,'children'),'ylim',[min_ylim max_ylim]);
        SetFigure(15);
        
    end

%%%%%%%%%%%%%%%%%%%  6. SVM Parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ SVM Training ------
select_for_SVM = select_bottom_line;
j_for_SVM = 1;
min_reps_training = 15; % For screening cells. (Anne used 20)
training_reps = 25; % I artificially added some trials by bootstrapping if reps < training_reps.
% This can include all cells (as Gu suggested) while not change the performance too much. @HH20150521
SVM_training_epoch = [700 1000]; % Anne used 500~700 ms
n_bootstrap_weights = 1000; % Anne used 1000
n_bootstrap_threshold = 50; % Anne used 25

% ------ SVM Testing ------
SVM_testing_span = 100; % Anne used 100 ms
SVM_testing_step = 50;
% n_bootstrap_test = 300; % Anne used 1000
n_bootstrap_test = 1000; % Anne used 1000


% ------ Initialization -----
pseudo_trial_pool_decoder = [];  pseudo_trial_pool_testing = [];
SVM_testing_tcenters = min(rate_ts{j_for_SVM}) + SVM_testing_span/2 : SVM_testing_step : max(rate_ts{j_for_SVM}) - SVM_testing_span/2;
answers_choice = []; answers_modality = []; reps_actual_training = [];
weights_svm_choice_allbootstrap = []; weights_svm_modality_allbootstrap = []; angles = [];
weights_svm_choice_mean = []; weights_svm_modality_mean = [];
thres_choice = []; thres_modality = []; select_for_SVM_actual = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dPCA_data;

    function f5p4(debug)    % dPCA is finally added by ZZ @20220310
        if debug;  dbstack;  keyboard;  end
%                 selected_dPCA_data = dPCA_data(select_tcells);
        selected_dPCA_data = dPCA_data(select_cpref_mpref);
%         for st = 1:3
%             %         dPCA_ZZ(group_result, CORRECTNWRONG_ANGLE,st);
%             dPCA_ZZ(selected_dPCA_data, CORRECT_ANGLE,st);     % For each condition, sorted by headings and choice
%         end
        
        dPCA_ZZ(selected_dPCA_data, ALL_CorrectCHOICE,[]);        % For all conditions, sorted by modality and choice
        
%         % Gouki's suggestion: trajectory of three modalities in the same dPCA subspace
%         dPCA_ZZ(selected_dPCA_data, CORRECT_ANGLE, []); 
        
    end

    function f5p5(debug)    % TDR added by ZZ @20231205
        if debug;  dbstack;  keyboard;  end
        
        selected_tdr_data = dPCA_data(select_cpref_mpref);
        
        TDR_ZZ(selected_tdr_data); 
        
    end


    function f6p0(debug)      % SVM 1: Training @HH20150429
        if debug;  dbstack;  keyboard;  end
        
        
        % Override
        %         min_reps = 15; % Anne used 20
        %         training_reps = 25; % I artificially added some trials by bootstrapping if reps < training_reps.
        %                             % This can include all cells (as Gu suggested) while not change the performance too much. @HH20150521
        
        %         SVM_training_epoch = [0 1600];    % This make it comparable with PCA results
        %         SVM_training_epoch = [1500 1700];    % This make it comparable with PCA results
        
        find_for_SVM = find(select_for_SVM);
        SVM_training_epoch_ind = min(SVM_training_epoch) <= rate_ts{j_for_SVM} & rate_ts{j_for_SVM} <= max(SVM_training_epoch);
        
        % Regroup pseudo-trials pool
        pseudo_trial_pool_decoder = cell(length(find_for_SVM),4);
        pseudo_trial_pool_testing = cell(length(find_for_SVM),4,length(SVM_testing_tcenters));
        
        for ii = 1:length(find_for_SVM)
            raw_this = group_result(find_for_SVM(ii)).mat_raw_PSTH.PSTH{j_for_SVM,ALL_CorrectCHOICE,1}.raw([1 2 3 4]); % Vest & vis
            %             raw_this = group_result(find_for_SVM(ii)).mat_raw_PSTH.PSTH{j_for_SVM,ALL_CHOICE,1}.raw([1 2 5 6]); % Vest & Comb
            %                         raw_this = group_result(find_for_SVM(ii)).mat_raw_PSTH.PSTH{j_for_SVM,ALL_CHOICE,1}.raw([3 4 5 6]); % Vis & Comb
            
            % By default, SVM for choice is related to pref and null
            % This make it comparable with PCA results and reasonable for TDR
            
            % Flip to contralateral vs ipsilateral
            %             if ~group_result(find_for_SVM(ii)).if_contralateral
            %                 raw_this = raw_this([2 1 4 3]);  % +/-,+/-
            %             end
            
            % Flip to right vs left
            %             if group_result(find_for_SVM(ii)).PREF_PSTH == 1
            %                 raw_this = raw_this([2 1 4 3]);
            %             end
            
            % Training pool (one certain time window)
            means_this = cellfun(@(x)mean(x(:,SVM_training_epoch_ind),2),raw_this,'uniformOutput',false);
            [pseudo_trial_pool_decoder{ii,:}] = means_this{:};
            
            % Testing pool (shifting time windows)
            for ttt = 1:length(SVM_testing_tcenters)
                test_epoch_ind = SVM_testing_tcenters(ttt) - SVM_testing_span/2 <= rate_ts{j_for_SVM} & rate_ts{j_for_SVM} <= SVM_testing_tcenters(ttt) + SVM_testing_span/2;
                means_this = cellfun(@(x)mean(x(:,test_epoch_ind),2),raw_this,'uniformOutput',false);
                
                [pseudo_trial_pool_testing{ii,:,ttt}] = means_this{:};
            end
            
        end
        
        % Find cells who cross the minimal repetitions for each condition (Anne)
        lengths = cellfun(@(x)length(x),pseudo_trial_pool_decoder);
        select_for_SVM_actual = all(lengths >= min_reps_training,2);
        pseudo_trial_pool_decoder = pseudo_trial_pool_decoder(select_for_SVM_actual,:);
        pseudo_trial_pool_testing = pseudo_trial_pool_testing(select_for_SVM_actual,:,:);
        
        lengths = lengths(select_for_SVM_actual,:);
        
        % I artificially added some trials by bootstrapping if reps < training_reps.
        % This can include all cells (as Gu suggested) while not change the performance too much. @HH20150521
        
        who_needs_add_some_trials = find(lengths < training_reps);
        
        for aa = 1:length(who_needs_add_some_trials)
            this_who = who_needs_add_some_trials(aa);
            
            % Select trials
            real_reps = length(pseudo_trial_pool_decoder{this_who});
            n_to_add = training_reps - real_reps;
            trials_to_add = randperm(real_reps,n_to_add);
            
            % Add trials
            pseudo_trial_pool_decoder{this_who} = [pseudo_trial_pool_decoder{this_who}; pseudo_trial_pool_decoder{this_who}(trials_to_add)];
            [sub1,sub2] = ind2sub(size(pseudo_trial_pool_decoder),this_who);
            
            for ttt = 1:length(SVM_testing_tcenters)
                pseudo_trial_pool_testing{sub1,sub2,ttt} = [pseudo_trial_pool_testing{sub1,sub2,ttt}; pseudo_trial_pool_testing{sub1,sub2,ttt}(trials_to_add)];
            end
        end
        
        % Recalculate lengths
        lengths = cellfun(@(x)length(x),pseudo_trial_pool_decoder);
        reps_actual_training = min(lengths);
        
        % -------- Teacher signals --------
        % For choices, 1 = Contralateral, 0 = Ipsilateral
        answers_choice = [ones(1,reps_actual_training(1)) zeros(1,reps_actual_training(2)) ones(1,reps_actual_training(3)) zeros(1,reps_actual_training(4))]';
        % For modality, 0 = Vestibular, 1 = Visual
        answers_modality = [zeros(1,reps_actual_training(1)) zeros(1,reps_actual_training(2)) ones(1,reps_actual_training(3)) ones(1,reps_actual_training(4))]';
        
        %% =============  SVM Bootstrapping ============
        % Parallel computing
        if ver_num < 2014
            if matlabpool('size') == 0 ;  matlabpool;  end
        else
            if isempty(gcp('nocreate')); parpool; end
        end
        
        training_set = cell(n_bootstrap_weights);
        tic; disp('Starting SVM training...');
        
        parfor_progress(n_bootstrap_weights);
        
        parfor nn = 1:n_bootstrap_weights
            % Randomly select trials without replacement from the pseudo-trials pool
            pseudo_trial_pool_perm = cellfun(@(x)x(randperm(size(x,1)),1),pseudo_trial_pool_decoder,'uniform',false);
            
            for conds = 1:4
                training_set{nn}{conds} = cell2mat(cellfun(@(x)x(1:reps_actual_training(conds),:),pseudo_trial_pool_perm(:,conds),'uniform',false)');
            end
            
            firing_for_this_training = cell2mat(training_set{nn}');
            
            % Train two SVM classifiers (choice and modality)
            svm_training_choice(nn) = svmtrain(firing_for_this_training,answers_choice,'box',1e-5,'tol',1e-7,'autoscale',0);
            svm_training_modality(nn) = svmtrain(firing_for_this_training,answers_modality,'box',1e-5,'tol',1e-7,'autoscale',0);
            
            parfor_progress;
        end
        parfor_progress(0);
        toc; disp('Done!');
        
        %% Averaged weights (bagging)
        %         weights_svm_choice_allbootstrap = cell2mat(arrayfun(@(x)-x.SupportVectors'*x.Alpha,svm_training_choice,'uniform',0));
        %         weights_svm_modality_allbootstrap = cell2mat(arrayfun(@(x)-x.SupportVectors'*x.Alpha,svm_training_modality,'uniform',0));
        
        % By default, 'autoscale' = true, so here I should rescale back the weight! HH20170613
        weights_svm_choice_allbootstrap = cell2mat(arrayfun(@(x)-x.SupportVectors'*x.Alpha,svm_training_choice,'uniform',0));
        weights_svm_modality_allbootstrap = cell2mat(arrayfun(@(x)-x.SupportVectors'*x.Alpha,svm_training_modality,'uniform',0));
        
        % Orthogonal test
        angles = acos(sum(weights_svm_choice_allbootstrap.*weights_svm_modality_allbootstrap,1)./...
            sqrt(sum(weights_svm_choice_allbootstrap.^2,1)./ sum(weights_svm_modality_allbootstrap.^2,1)))/pi*180;
        
        weights_svm_choice_mean = mean(weights_svm_choice_allbootstrap,2);
        weights_svm_choice_mean = weights_svm_choice_mean / max(abs(weights_svm_choice_mean));
        % svm_weights_choice_sem = std(svm_weights_choice,[],2)/sqrt(n_bootstrap);
        
        weights_svm_modality_mean = mean(weights_svm_modality_allbootstrap,2);
        weights_svm_modality_mean = weights_svm_modality_mean / max(abs(weights_svm_modality_mean));
        % svm_weights_modality_sem = std(svm_weights_modality,[],2)/sqrt(n_bootstrap);
        
        %% Find proper threshold for SVM
        training_set = cell(n_bootstrap_threshold);
        bestZ_choice = nan(n_bootstrap_threshold,1);
        bestZ_modality = nan(n_bootstrap_threshold,1);
        tic; disp('Find proper threshold for SVM decoders...');
        
        parfor nn = 1:n_bootstrap_threshold
            % Randomly select trials without replacement from the pseudo-trials pool
            pseudo_trial_pool_perm = cellfun(@(x)x(randperm(size(x,1)),1),pseudo_trial_pool_decoder,'uniform',false);
            
            for conds = 1:4
                training_set{nn}{conds} = cell2mat(cellfun(@(x)x(1:reps_actual_training(conds),:),pseudo_trial_pool_perm(:,conds),'uniform',false)');
            end
            
            proj_choice_on_choice = cell2mat(training_set{nn}') * weights_svm_choice_mean;
            proj_modality_on_modality = cell2mat(training_set{nn}') * weights_svm_modality_mean;
            
            % Find the best threshold (the most northwestern point on ROC curve)
            [~,bestZ_choice(nn)] = rocN(proj_choice_on_choice(logical(answers_choice)),proj_choice_on_choice(~logical(answers_choice)));
            [~,bestZ_modality(nn)] = rocN(proj_modality_on_modality(logical(answers_modality)),proj_modality_on_modality(~logical(answers_modality)));
            
        end
        toc; disp('Done!');
        
        thres_choice = mean(bestZ_choice);
        thres_modality = mean(bestZ_modality);
        
        % ================ SVM Training End ===============
        
    end

    function f6p1(debug)      % SVM 2: Plotting SVM weights
        if debug;  dbstack;  keyboard;  end
        
        
        if isempty(thres_choice)
            f6p0(0);
        end
        
        % ------ Plotting -----
        set(figure(611),'pos',[89 26 900 800]); clf
        
        find_for_SVM = find(select_for_SVM);
        find_for_SVM_actual_of_all = find_for_SVM(select_for_SVM_actual);
        who_are_tcells = select_tcells(find_for_SVM_actual_of_all);
        
        % Weights sorted by choice decoder added by ZZ 20210514
        % Weights for choice decoder
        subplot(3,2,1);
        [~,sortt] = sort(abs(weights_svm_choice_mean),'descend');
        set(bar(weights_svm_choice_mean(sortt),0.5),'facecolor','k','edgecolor','none'); hold on;
        set(bar(find(who_are_tcells(sortt)),...
            weights_svm_choice_mean(sortt(who_are_tcells(sortt))),0.5),'facecolor',condition_colors(2,:),'edgecolor','none');
        
        xlim([-2 sum(select_for_SVM_actual)+1]); ylim([-1 1]); title('Choice weights');
        set(gca,'xtick',[1 sum(select_for_SVM_actual)]);
        % errorbar(svm_weights_choice_mean(sortt),svm_weights_choice_sem(sortt),'.');
        
        subplot(3,2,2);
        set(bar(weights_svm_modality_mean(sortt), 0.5), 'facecolor', 'k', 'edgecolor', 'none'); hold on;
        set(bar(find(who_are_tcells(sortt)),...
            weights_svm_modality_mean(sortt(who_are_tcells(sortt))), 0.5), 'facecolor', condition_colors(2,:), 'edgecolor', 'none');
        xlim([-2 sum(select_for_SVM_actual)+1]);  ylim([-1 1]); title('Modality weights sorted by choice weights');
        set(gca,'xtick',[1 sum(select_for_SVM_actual)]);
        
        % Weights sorted by modality decoder
        % Weights for modality decoder
        subplot(3,2,3);
        [~,sortt] = sort(abs(weights_svm_modality_mean),'descend');
        set(bar(weights_svm_modality_mean(sortt),0.5),'facecolor','k','edgecolor','none'); hold on;
        set(bar(find(who_are_tcells(sortt)),...
            weights_svm_modality_mean(sortt(who_are_tcells(sortt))),0.5),'facecolor',condition_colors(2,:),'edgecolor','none');
        
        xlim([-2 sum(select_for_SVM_actual)+1]);  ylim([-1 1]); title('Modality weights');
        set(gca,'xtick',[1 sum(select_for_SVM_actual)]);
        % errorbar(svm_weights_modality_mean(sortt),svm_weights_modality_sem(sortt),'.');
        
        subplot(3,2,4);
        set(bar(weights_svm_choice_mean(sortt),0.5), 'facecolor', 'k', 'edgecolor', 'none'); hold on;
        set(bar(find(who_are_tcells(sortt)),...
            weights_svm_choice_mean(sortt(who_are_tcells(sortt))), 0.5), 'facecolor', condition_colors(2,:), 'edgecolor', 'none');
        
        xlim([-2 sum(select_for_SVM_actual)+1]); ylim([-1 1]); title('Choice weights sorted by madality weights');
        set(gca,'xtick',[1 sum(select_for_SVM_actual)]);
        
        
        
        % Correlation
        monkeys = xls_num{1}(:,header.Monkey);  % Temporarily here. HH20150914
        monkey1 = monkeys == 15;
        
        %         monkey1 = monkeys == 5;
        %         monkey2 = monkeys == 10;
        
        ax = subplot(3,2,5);
        h = LinearCorrelation({
            (weights_svm_modality_mean(monkey1(select_for_SVM)));
            %             (weights_svm_modality_mean(monkey2(select_for_SVM)))
            },...
            {
            (weights_svm_choice_mean(monkey1(select_for_SVM))) ;
            %             (weights_svm_choice_mean(monkey2(select_for_SVM)))
            },...
            'CombinedIndex',[],...
            'Xlabel','Weights in modality decoder','Ylabel','Weights in choice decoder',...
            'FaceColors',{'k'},'Markers',{'o'},...
            'LineStyles',{'k-'},'MarkerSize',6,...
            'axes',ax,'XHist',0,'YHist',0,...
            'SameScale',1,'Method','Pearson','FittingMethod',2);
        %         delete([h.diag h.group(1:2).line]);
        legend off;
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');
        %         text(min(xlim)+0.2,min(ylim)+0.2,sprintf('r^2 = %g, p = %g',h.group(3).r_square,h.group(3).p),'fontsize',11);
        
        % Annotate tcells
        h_line = plot(weights_svm_modality_mean(who_are_tcells),...
            weights_svm_choice_mean(who_are_tcells),...
            '+','markersize',8,'color',condition_colors(2,:),'linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        select_for_SVM_actual_of_all = zeros(N,1);
        select_for_SVM_actual_of_all(find_for_SVM_actual_of_all) = 1;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h.group(1).dots, select_for_SVM_actual_of_all});
        
        % Angles between weights
        subplot(3,2,6);
        [n,x] = hist(angles,20); bar(x,n,'facecolor','k','edgecolor','k'); hold on;
        plot([median(angles) median(angles)],ylim,'k-');
        plot([90 90],ylim,'k--'); axis tight;
        %         xlim([0 180]); set(gca,'xtick',0:90:180);
        title(sprintf('Median angle = %g',median(angles)));
        
        SetFigure(15);
        
        % ---------  Memsac DDI and weights for choice decoder
        set(figure(612),'pos',[89 26 761 606]); clf
        
        MemSac_DDI_phase = [2:4];
        
        h = LinearCorrelation({
            MemSac_indicator(select_for_SVM);
            },...
            {
            (weights_svm_choice_mean);
            },...
            'CombinedIndex',[],...
            'Xlabel',MemSac_indicator_txt,'Ylabel','Weight in choice decoder',...
            'FaceColors',{'k'},'Markers',{'o'},'Markers',{'o'},...
            'LineStyles',{'k-'},'MarkerSize',12,...
            'figN',612,'XHist',20,'YHist',20,...
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Pearson','FittingMethod',2);
        
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        
        % Annotate tcells
        plot(MemSac_indicator(select_tcells),weights_svm_choice_mean(select_tcells(select_for_SVM),1),...
            '+','markersize',16,'color',condition_colors(2,:),'linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot(MemSac_indicator(select_for_SVM),(weights_svm_choice_mean(:,1)),'visible','off'); hold on;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_for_SVM});
        
        
        % ---------  Weights for eigen-neuron 1 and weights for SVM choice decoder
        if isempty(weights_PCA_B_PC)
            f5p0(0);
        end
        
        set(figure(613),'pos',[89 26 761 606]); clf
        
        h = LinearCorrelation({
            Modality_pref_all(1,:,1);
            % weights_PCA_B_PC(:,1);
            },...
            {
            % (weights_svm_choice_mean);
            (weights_svm_modality_mean);
            },...
            'CombinedIndex',[],... % 'Xlabel','Weight for eigen-neuron 1','Ylabel','Weight in choice decoder',...
            'Xlabel','Modality pref (vis - vest)','Ylabel','Weight in modality decoder',...
            'FaceColors',{'k'},'Markers',{'o'},'Markers',{'o'},...
            'LineStyles',{'k-'},'MarkerSize',12,...
            'figN',613,'XHist',20,'YHist',20,...
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Spearman','FittingMethod',2);
        
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        
        %         % Annotate tcells
        %         plot(weights_PCA_B_PC(select_tcells(select_for_SVM),1),weights_svm_choice_mean(select_tcells(select_for_SVM),1),...
        %             '+','markersize',16,'color',colors(2,:),'linew',2);
        %
        %         % Show individual cell selected from the figure. HH20150424
        %         set([gca h.group.dots],'ButtonDownFcn',{@Show_individual_cell, h.group.dots, select_for_SVM});
        
        % ------------ Weights for SVM choice decoder VS choice divergence and multisensory enhancement 20170609 ------------------
        %         Choice_pref_all_comb = abs(Choice_pref_all(3,:,3)); % stim-on to stim-off
        %         Choice_pref_all_vest = abs(Choice_pref_all(1,:,3));
        %         Choice_pref_all_vis = abs(Choice_pref_all(2,:,3));
        Choice_pref_all_comb = abs(Choice_pref_all(3,:,1)); % stim-on to stim-off
        Choice_pref_all_vest = abs(Choice_pref_all(1,:,1));
        Choice_pref_all_vis = abs(Choice_pref_all(2,:,1));
        
        figure(1234); clf;
        hs = tight_subplot(2,2,0.1,0.1,0.1);
        
        % -- SVM weights vs Div_vest --
        xxs = {Choice_pref_all_vest(:), 'Divergence vest';
            Choice_pref_all_vis(:),  'Divergence vis';
            Choice_pref_all_comb(:), 'Divergence comb';
            % Choice_pref_all_comb(:)./(Choice_pref_all_vest(:)+Choice_pref_all_vis(:)), '(comb-vest)+(comb-vis)'};
            % 2*Choice_pref_all_comb(:)-(Choice_pref_all_vest(:)+Choice_pref_all_vis(:)), '(comb-vest)+(comb-vis)'};
            Choice_pref_all_comb(:)-(Choice_pref_all_vest(:)), '(comb-vest)+(comb-vis)'};
        
        for xxxx = 1:length(xxs)
            xx = xxs{xxxx,1};
            yy = (weights_svm_choice_mean);
            hl = LinearCorrelation(xx,yy,'Axes',hs(xxxx),'MethodOfCorr','Spearman','FittingMethod',2);
            hold on;
            delete([hl.leg]);
            
            xlabel(xxs{xxxx,2}); ylabel('abs(svm weight)');
            text(min(xlim),min(ylim)+range(ylim)*0.9,...
                sprintf('r^2=%g\n p=%g\n k = %g\\pm%g',hl.group.r_square,hl.group.p,hl.group.para(1),hl.group.paraSE(1)))
            
            % Annotate tcells
            plot(xx(select_tcells),yy(select_tcells(select_for_SVM),1),...
                '+','markersize',16,'color',condition_colors(2,:),'linew',2);
            
        end
        SetFigure(20);
        
        
    end
    function f6p2(debug)      % SVM 3: Testing SVM
        if debug;  dbstack;  keyboard;  end
        
        
        if isempty(thres_choice)
            f6p0(0);
        end
        
        corr_rate_choice_by_choice = nan(n_bootstrap_test,length(SVM_testing_tcenters));
        corr_rate_modality_by_modality = corr_rate_choice_by_choice;
        corr_rate_choice_by_modality = corr_rate_choice_by_choice;
        corr_rate_modality_by_choice = corr_rate_choice_by_choice;
        
        tic; disp('Testing SVM decoders...');
        
        progressbar('Testing SVM');
        
        for nn = 1:n_bootstrap_test  % For each bootstrapping
            % Generate rand perm indices
            pseudo_trial_pool_perm_ind = cellfun(@(x)randperm(size(x,1)),pseudo_trial_pool_decoder,'uniform',false);
            
            for ttt = 1:length(SVM_testing_tcenters)
                
                for conds = 1:4
                    for ii = 1:size(pseudo_trial_pool_perm_ind,1)
                        testing_set{nn}{conds}(:,ii) = pseudo_trial_pool_testing{ii,conds,ttt}(pseudo_trial_pool_perm_ind{ii,conds}(1:reps_actual_training(conds)));
                    end
                end
                
                % Projection on SVM weights
                proj_choice_on_choice = cell2mat(testing_set{nn}') * weights_svm_choice_mean;
                proj_modality_on_modality = cell2mat(testing_set{nn}') * weights_svm_modality_mean;
                
                % Classification by thresholds
                class_test_choice = proj_choice_on_choice >= thres_choice;
                class_test_modality = proj_modality_on_modality >= thres_modality;
                
                % Correct rates
                corr_rate_choice_by_choice(nn,ttt) = sum(answers_choice == class_test_choice)/length(answers_choice);
                corr_rate_modality_by_modality(nn,ttt) = sum(answers_modality == class_test_modality)/length(answers_modality);
                corr_rate_choice_by_modality(nn,ttt) = sum(answers_choice == class_test_modality)/length(answers_choice);
                corr_rate_modality_by_choice(nn,ttt) = sum(answers_modality == class_test_choice)/length(answers_modality);
                
            end
            
            progressbar(nn/n_bootstrap_test);
        end
        
        toc; disp('Done!');
        
        %% --------- Plotting ---------
        
        set(figure(621),'pos',[114 223 930 353]); clf
        subplot(1,2,1);
        plot(SVM_testing_tcenters,nanmean(corr_rate_choice_by_choice), 'color', 'm', 'LineWidth', 2); hold on;
        plot(SVM_testing_tcenters,nanmean(corr_rate_modality_by_choice), 'color', [.87 .49 0],'LineWidth', 2);
        legend('Decodes choice','Decodes modality','location','best', 'AutoUpdate', 'off');
        
        shadedErrorBar(SVM_testing_tcenters,nanmean(corr_rate_choice_by_choice),...
            nanstd(corr_rate_choice_by_choice),'lineprops',{'color','m','linew',2});
        shadedErrorBar(SVM_testing_tcenters,nanmean(corr_rate_modality_by_choice),...
            nanstd(corr_rate_modality_by_choice),'lineprops',{'color',[.87 .49 0],'linew',2});
        
        axis tight;  ylim([0.2 1.05]);
        plot(xlim,[0.5 0.5],'k--'); title('Choice decoder');
        
        % Time markers
        for tt = 1:2
            plot([1 1] * time_markers{j_for_SVM}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j_for_SVM}{tt},'linew',1);
        end
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{j_for_SVM}(1), .5+ Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
        
        subplot(1,2,2);
        shadedErrorBar(SVM_testing_tcenters,nanmean(corr_rate_modality_by_modality),...
            nanstd(corr_rate_modality_by_modality),'lineprops',{'color',[.87 .49 0],'linew',2}); hold on;
        shadedErrorBar(SVM_testing_tcenters,nanmean(corr_rate_choice_by_modality),...
            nanstd(corr_rate_choice_by_modality),'lineprops',{'color','m','linew',2});
        
        axis tight; ylim([0.2 1.05]);
        plot(xlim,[0.5 0.5],'k--'); title('Modality decoder');
        
        % Time markers
        for tt = 1:2
            plot([1 1] * time_markers{j_for_SVM}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j_for_SVM}{tt},'linew',1);
        end
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{j_for_SVM}(1),.5+ Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
        
        SetFigure(15);
        
        %%
        %         keyboard
    end
    function f6p3(debug)      % SVM 4: Mean PSTH projectd on weight vector (SVM-weighted sum)
        if debug;  dbstack;  keyboard;  end
        
        
        if isempty(thres_choice)
            f6p0(0);
        end
        
        find_for_SVM = find(select_for_SVM);
        find_for_SVM_actual_of_all = find_for_SVM(select_for_SVM_actual);
        select_for_SVM_actual_of_all = false(N,1);
        select_for_SVM_actual_of_all(find_for_SVM_actual_of_all) = true;
        
        % projs = {weights_svm_choice_mean,weights_svm_modality_mean};
        projs = {weights_svm_choice_allbootstrap,weights_svm_modality_allbootstrap}; % With bootstrap
        selects = {select_for_SVM_actual_of_all,select_for_SVM_actual_of_all};
        
        h = Weighted_sum_PSTH(projs,{'Weighted by svm\_choice\_mean','Weighted by svm\_modality\_mean'},selects);
        set(h.fig,'name','Projected on SVM decoders (Correct only, all choices)');
        
    end


% Added by ZZ @ 20220228
% Lasso Decoder training
% Refer to Kiani et al., 2014, Cur. Biolog.
% Logit(P_right(t)) = Beta0(t) + Beta(t) .* FR(t) ;
% Beta is a vector with n*1 dimension, n is the number of neuron
% Changed @ 20230705

% % ------------------  Pre-determination  ---------------
% select_for_Lasso = select_bottom_line;
% select_for_lasso_decoder = [];
% j_for_Lasso = [1;2];
% min_reps_training = 40;
% min_reps_decoder = 50;
% % training_reps = 25; % I artificially added some trials by bootstrapping if reps < training_reps.
% Lasso_sliding_window = 100;  % ms, from Kiani et al., 2014
% Lasso_step_size = 20the sliding ;
% Lasso_bootstrap = 100;
% 
% for j = 1:length(j_for_Lasso)
%     Lasso_tcenters{j} = min(rate_ts{j_for_Lasso(j)}) + Lasso_sliding_window/2 : Lasso_step_size : max(rate_ts{j_for_Lasso(j)}) - Lasso_sliding_window/2;
%     pseudo_trial_lasso_decoder{j} = cell(length(select_for_Lasso),6,length(Lasso_tcenters{j}));
%     
% end
% 
% 
% % Initiation
% pseudo_trial_pool_decoder = cell(length(j_for_Lasso),1);
% pseudo_trial_pool_perm = cell(length(j_for_Lasso),Lasso_bootstrap);
% training_set = cell(length(j_for_Lasso),Lasso_bootstrap);
% B = cell(length(Lasso_tcenters{j}),2,Lasso_bootstrap,size(pseudo_trial_lasso_decoder{j},2)/2);
% FitInfo = cell(length(Lasso_tcenters{j}),2,Lasso_bootstrap,size(pseudo_trial_lasso_decoder{j},2)/2);
% reps_actual_decoder =[];
model = []; 
decoder_data = []; 
accuracy_training = [];  % only for svm
accuracy_test = []; 
SVM = 1;   LASSO_A = 2;  % SVM for comparison with Lasso

% Data (Training and testing data) Preparing
    function f11p0(debug)
        if debug;  dbstack;  keyboard;  end
        
        bootstrap_n = 1000;
        % Select data for decoders
        select_lasso_data = dPCA_data(select_cpref_mpref);
        
        if isempty(decoder_data)
            if exist('decoder_data.mat', 'file')
                load('decoder_data.mat')
            else
            % == Preparing manually concatenated population data for Training and testing ==
            decoder_data = decoder_training_prepare(select_lasso_data,All_Choice,'bootstrapN',bootstrap_n);
            end
        end
        
    end
        
% Model training 
    function f11p1(debug)
        if debug;  dbstack;  keyboard;  end
        
        if isempty(decoder_data)
            f11p0(0);
        end
        
        if isempty(model)
            if exist('model.mat', 'file')
                load('model.mat')
            else
                % == Training ==
                %             for type = SVM:LASSO_A
                progressbar( 'Lasso training' );
                for type = LASSO_A
                    for st = 1:3
                        model{type, st} = decoder_train(decoder_data.FR,decoder_data.training_set, decoder_data.teacher_signal,st,type);
                        
                        progressbar(st/3);
                    end
                end
            end
        end
    end
        

% Visualing decoders 
    function f11p2(debug)
        if debug;  dbstack;  keyboard;  end
        
        if isempty(model)
            f11p1(0); 
        end
        
        % == Visualize the weight of decoders 
        progressbar('Plot Decoder Weight')
        for type = LASSO_A
%         for type = SVM : LASSO_A
            for st = 1:3
                plot_decoder_correlation(model{type,st}, type, st, decoder_data.t_centers);
                
                progressbar(((type-1)*3+st) / 6);
            end
        end
        
    end


% Across-modality decoding, only for Lasso
    function f11p3(debug)
        if debug;  dbstack;  keyboard;  end
        
%         save('decoder_prepare_data.mat', 'decoder_data', '-v7.3');
%         save('trained_decoder.mat', 'model', '-v7.3');
%         save('decoder_tcenters.mat', 'decoder_tcenters', '-v7.3'); 
        
        if isempty(model)
            f11p1(0);
        end

        type = LASSO_A;
        progressbar('Cross-modal decoding')
        
        for  st_train = 1:3
            for st_test = 1:3
                accuracy_training{st_train, st_test} = decoder_test(model{type,st_train},decoder_data.FR, decoder_data.training_set,st_test,type);
                
                progressbar(((st_train-1)*3 +st_test)/9);
            end
        end
        
    end


    function f11p4(debug)
        if debug;  dbstack;  keyboard;  end
        
        if isempty(accuracy_training)
            f11p3(0);
        end

        type = LASSO_A;
        
        for st_train = 1:3
            for st_test = 1:3
                plot_temporal_decoder_accuracy(accuracy_training{st_train,st_test}, type, st_test, decoder_data.t_centers);
                
            end
        end
    end


% Within modality decoding
    function f11p5(debug)
        if debug;  dbstack;  keyboard;  end
        
        if isempty(model)
            f11p1(0);
        end

        % == Testing == 
        for type = LASSO_A
%         for type = SVM:LASSO_A
            progressbar('Decoder Testing')
            for st = 1:3
%                 accuracy{type, st} = decoder_test(model{type,st},decoder_data.training_set,st,type);
                accuracy_test{type, st} = decoder_test(model{type,st},decoder_data.FR, decoder_data.testing_set,st,type);
                progressbar(st/3);
            end
        end
    end
        

    function f11p6(debug)
        if debug;  dbstack;  keyboard;  end
        
        if isempty(accuracy_test)
            f11p5(0);
        end
        
        % Visualize performance of decoders
        progressbar('Plot Performance of Decoder')
        for type = LASSO_A
%         for type = SVM:LASSO_A
            for st = 1:3
%                 plot_temporal_decoder_accuracy(accuracy{type,st}, type, st, decoder_tcenters);
                plot_temporal_decoder_accuracy(accuracy_test{type,st}, type, st,decoder_data.t_centers);
                
                progressbar(((type-1)*3+st) / 6);
            end
        end
    end

fisherSimpleGu = [];
findForFisher = select_cpref_mpref + select_tcells;    % Changed by ZZ for the exclusion of 33 outlier cells 
findForFisher(findForFisher==0) = [];
findForFisher = find(findForFisher==2);

%% ====== Calculate Fisher information of heading ======== HH20180619
% 1. Simplest method like Gu 2010: Sum over cells (slope / mean)
    function f6p5p1(debug)
        if debug;  dbstack;  keyboard;  end
        %%
        
        j = 1;
        
        fisherSize = size(group_result(1).mat_raw_PSTH.fisherSimpleGu);
        fisherSimpleGu = nan(length(findForFisher),fisherSize(1),fisherSize(2),fisherSize(3),fisherSize(4));  % Last one: varPoisson/varReal
        
        for ii = 1:length(findForFisher)
            
            thisCell = findForFisher(ii);
            fisherSimpleGu (ii,:,:,:,:) = group_result(thisCell).mat_raw_PSTH.fisherSimpleGu;
            
        end
        
        % -------- Plotting FI(t) ---------
        figure(4232); clf
        set(gcf,'uni','norm','pos',[0.442       0.166       0.557       0.739]);
        %         plotRange = 1:find(CP_ts{1}>=1500,1);
        plotRange = 1: size(CP_ts{1},2);
        
        bootN = 1000;
        
        rangeTitles = {'+/-8 degree','+/-2 degrees'};
        titles = {'Poisson assump.','Real variance'};
        
        for ranges = 1:2
            for varMethod = 1:2  %  1: Gu's Poisson assumption   2: Real variance (10x larger than Poisson because of varCE)
                h = subplot(2,2,varMethod + (ranges-1)*2);
                
                % Get sum of Fisher and std by bootstrap (Gu 2010)
                %                 sumBoots = bootstrp(bootN,@(x)sum(x,1),squeeze(fisherSimpleGu(:,plotRange,:,ranges,varMethod)));
                sumBoots = bootstrp(bootN,@(x)nansum(x,1),squeeze(fisherSimpleGu(:,plotRange,:,ranges,varMethod)));
                sumBoots = reshape(sumBoots,bootN,[],3);
                
                SeriesComparison(sumBoots, CP_ts{1}(plotRange),...
                    'SEM',0,'Errorbar',2,'Axes',h,'Transparent',transparent);
                
                sumFisherMean = squeeze(mean(sumBoots,1));
                sumFisherVestPlusVIs = sumFisherMean(:,1) + sumFisherMean(:,2);
                plot(CP_ts{1}(plotRange),sumFisherVestPlusVIs,'m','linew',2);
                
                % Gaussian vel
                axis tight; plot([0 0],ylim,'k--'); plot([1500 1500], ylim,'k--')
                xlim([-100 1600]); ylim([0 max(ylim)*1.1])
                plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2.5,'color',[0.6 0.6 0.6]);
                legend off;
                title(titles{varMethod})
                ylabel(rangeTitles{ranges});
            end
        end
        
        SetFigure(15)
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FIScatterTimeRange = 500 < CP_ts{1} & CP_ts{1} <= 1500;
FIScatterHeadingRange = 2; % +/-8 degrees
FIScatterVariance = 2; % Real variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1.1 Plot scatter FI of different modalities
    function f6p5p1p1(debug)
        if debug;  dbstack;  keyboard;  end
        
        if isempty(fisherSimpleGu), f6p5p1(0); end
        
        % -------- Plotting cell-by-cell FI (Yong Gu wants this) ---------
        
        % Which info?
        FI_all_temp = squeeze(mean(fisherSimpleGu(:,FIScatterTimeRange,:,FIScatterHeadingRange,FIScatterVariance),2))';
        FI_temp_comb_minus_vest = FI_all_temp(3,:)-FI_all_temp(1,:);
        FI_temp_comb_minus_vis = FI_all_temp(3,:)-FI_all_temp(2,:);
        FI_temp_vest_plus_vis = FI_all_temp(1,:) + FI_all_temp(2,:);
        
        %%  1. ====== vest and visual  ========
        %         tt = 3;
        tt = 1;
        
        
        monkeys = xls_num{1}(:,header.Monkey);  % HH20160613
        monkey1 = monkeys == 15; monkey1 = monkey1(findForFisher)';
        monkey2 = monkeys == 13; monkey2 = monkey2(findForFisher)';
        
        set(figure(figN),'name','FI (visual) vs. FI (vest)','pos',[17 514 1151 449]);
        
        cpref_sig_1 = Choice_pref_p_value_all(1,findForFisher,tt) < 0.01;  % Use Choice pref p value as indicators
        cpref_sig_2 = Choice_pref_p_value_all(2,findForFisher,tt) < 0.01;
        
        
        h = LinearCorrelation({
            (FI_all_temp(2, monkey1 & ~cpref_sig_1 & ~cpref_sig_2)) ;
            (FI_all_temp(2, monkey2 & ~cpref_sig_1 & ~cpref_sig_2)) ;
            (FI_all_temp(2, monkey1 & xor(cpref_sig_1 , cpref_sig_2)));
            (FI_all_temp(2, monkey2 & xor(cpref_sig_1 , cpref_sig_2)));
            (FI_all_temp(2, monkey1 & cpref_sig_1 & cpref_sig_2));...
            (FI_all_temp(2, monkey2 & cpref_sig_1 & cpref_sig_2))
            },...
            {
            (FI_all_temp(1,monkey1 & ~cpref_sig_1 & ~cpref_sig_2)) ;
            (FI_all_temp(1,monkey2 & ~cpref_sig_1 & ~cpref_sig_2)) ;
            (FI_all_temp(1,monkey1 & xor(cpref_sig_1 , cpref_sig_2))) ;
            (FI_all_temp(1,monkey2 & xor(cpref_sig_1 , cpref_sig_2))) ;
            (FI_all_temp(1,monkey1 & cpref_sig_1 & cpref_sig_2)) ;...
            (FI_all_temp(1,monkey2 & cpref_sig_1 & cpref_sig_2))
            },...
            'CombinedIndex',[20 21 40 42],'PlotCombinedOnly', 1,...
            'Ylabel','Vestibular FI','Xlabel','Visual FI',...
            'FaceColors',{'none','none', [0.8 0.8 0.8],[0.8 0.8 0.8],'k','k'},'EdgeColors', {'k','k','k','k','k','k'}, 'Markers',{'o', '^'},...
            'LineStyles',{'k:','k:','k-'},'MarkerSize',marker_size,...
            'figN',figN,'XHist',20,'YHist',20,'Method','Pearson','FittingMethod',2, ...
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',1); figN = figN + 1;
        
        delete([h.group(1:2).line]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');       %  SetFigure(20);
        
        %         % Annotate tcells
        %         plot(FI_all_temp(2,select_tcells(select_tcells)),FI_all_temp(1,select_tcells(select_tcells)),...
        %             '+','markersize',tcell_cross_size,'color',colors(2,:),'linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot((FI_all_temp(2,:)),(FI_all_temp(1,:)),'visible','off'); hold on;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_tcells});
        
        
        %% 2. ====== Single vs Comb =======
        two_face_colors = fliplr({'none','none',[0.8 0.8 1],[0.8 0.8 1],condition_colors(1,:),condition_colors(1,:);
            'none','none',[1 0.8 0.8],[1 0.8 0.8],condition_colors(2,:),condition_colors(2,:)});
        
        for k = 1:2  % Plot it separately
            
            set(figure(figN),'name','FI (single) vs. FI (comb)','pos',[17 514 1151 449]);
            
            cpref_sig_k = Choice_pref_p_value_all(k,findForFisher,tt) < 0.05;
            %             cpref_sig_2 = Choice_pref_p_value_all(2,:,tt) < 0.05;
            cpref_sig_3 = Choice_pref_p_value_all(3,findForFisher,tt) < 0.05;
            
            
            h = LinearCorrelation({
                (FI_all_temp(k, monkey1 & cpref_sig_k & cpref_sig_3));
                (FI_all_temp(k, monkey2 & cpref_sig_k & cpref_sig_3));
                (FI_all_temp(k, monkey1 & xor(cpref_sig_k , cpref_sig_3)));
                (FI_all_temp(k, monkey2 & xor(cpref_sig_k , cpref_sig_3)));
                (FI_all_temp(k, monkey1 & ~cpref_sig_k & ~cpref_sig_3)) ;
                (FI_all_temp(k, monkey2 & ~cpref_sig_k & ~cpref_sig_3)) ;
                },...
                {
                (FI_all_temp(3,monkey1 & cpref_sig_k & cpref_sig_3)) ;
                (FI_all_temp(3,monkey2 & cpref_sig_k & cpref_sig_3));
                (FI_all_temp(3,monkey1 & xor(cpref_sig_k, cpref_sig_3))) ;
                (FI_all_temp(3,monkey2 & xor(cpref_sig_k, cpref_sig_3)));
                (FI_all_temp(3,monkey1 & ~cpref_sig_k & ~cpref_sig_3)) ;
                (FI_all_temp(3,monkey2 & ~cpref_sig_k & ~cpref_sig_3)) ;
                },...
                'CombinedIndex',[20 21 40 42],'PlotCombinedOnly', 1,...
                'Ylabel','Combined FI','Xlabel','Single FI',...
                'FaceColors',two_face_colors(k,:),'Markers',{'o','^'},...
                'EdgeColors', {'k','k','k','k','k','k'},...
                'LineStyles',{'k:','k:','k:','k:','k:','k:','k-'},'MarkerSize',marker_size,...
                'figN',figN,'XHist',20,'YHist',20,...
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',1,...
                'Method','Pearson','FittingMethod',2); figN = figN + 1;
            
            % delete([h.group(1:6).line]);
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');       %   SetFigure(20);
            axis square;
            
            % Annotate tcells
            %             plot((FI_all_temp(k,select_tcells(select_tcells))),(FI_all_temp(3,select_tcells(select_tcells))),...
            %                 '+','markersize',tcell_cross_size,'color','k','linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            h_line = plot((FI_all_temp(k,:)),(FI_all_temp(3,:)),'visible','off'); hold on;
            set([gca h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, [select_tcells],1});
            
            axis tight;
        end
        %             for i = 1:3
        %                 set(h.group(i).dots,'color',colors(k,:));
        %             end
        
        %% ===  3. (Comb - visual) VS (Comb - vest)  ===
        set(figure(figN),'name','FI(Comb - visual) VS FI(Comb - vest)','pos',[17 514 1151 449]);
        
        cellTypes = [group_result.Waveform_broad];
        cellTypes = cellTypes(findForFisher);
        cellTypes(:) = 0;
        
        
        h = LinearCorrelation({
            (FI_temp_comb_minus_vest(1,monkey1 & ~cellTypes)) ;
            (FI_temp_comb_minus_vest(1,monkey2 & ~cellTypes)) ;
            (FI_temp_comb_minus_vest(1,monkey1 & cellTypes)) ;
            (FI_temp_comb_minus_vest(1,monkey2 & cellTypes)) ;
            },...
            {
            (FI_temp_comb_minus_vis(1,monkey1 & ~cellTypes)) ;
            (FI_temp_comb_minus_vis(1,monkey2 & ~cellTypes)) ;
            (FI_temp_comb_minus_vis(1,monkey1 & cellTypes)) ;
            (FI_temp_comb_minus_vis(1,monkey2 & cellTypes)) ;
            },...
            'CombinedIndex',15,'PlotCombinedOnly', 1, ...
            'Xlabel','FI (Combined - Vest)','Ylabel','FI (Combined - Visual)',...
            'FaceColors',{'none','none','k','k'},'Markers',{'o','^'},'EdgeColors',{'k','k','k','k'},...
            'LineStyles',{'k:','k:','k:','k:','k-'},'MarkerSize',marker_size,...
            'figN',figN,... 'XHist',20,'YHist',20,'XHistStyle','stacked','YHistStyle','stacked',
            'SameScale',1,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        
        delete([h.group(1:2).line]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');         SetFigure(20);
        axis square;
        
        % Annotate tcells
        %         plot((FI_temp_comb_minus_vest(1,select_tcells(select_bottom_line))),...
        %             (FI_temp_comb_minus_vis(1,select_tcells(select_bottom_line))),...
        %             '+','markersize',tcell_cross_size,'color',colors(2,:),'linew',2);
        %             plot((Choice_pref_all_temp(2,select_tcells(select_bottom_line),tt)),(Choice_pref_all_temp(3,select_tcells(select_bottom_line),tt)),...
        %                 '+','markersize',tcell_cross_size,'color','k','linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot(FI_temp_comb_minus_vest(1,:),FI_temp_comb_minus_vis(1,:),'visible','off'); hold on;
        set([gca h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, [select_tcells],1});
        
        axis tight;
        
        
        %% ===  4. Comb VS (Vest + Vis)  ===
        set(figure(figN),'name','FI(Comb) VS FI(Vest + Vis)','pos',[17 514 1151 449]);
        
        cellTypes = [group_result.Waveform_broad];
        cellTypes(:) = 0;
        
        cellTypes = cellTypes(findForFisher);
        
        h = LinearCorrelation({
            (FI_temp_vest_plus_vis(1,monkey1 & ~cellTypes)) ;
            (FI_temp_vest_plus_vis(1,monkey2 & ~cellTypes)) ;
            (FI_temp_vest_plus_vis(1,monkey1 & cellTypes)) ;
            (FI_temp_vest_plus_vis(1,monkey2 & cellTypes)) ;
            },...
            {
            (FI_all_temp(3,monkey1 & ~cellTypes)) ;
            (FI_all_temp(3,monkey2 & ~cellTypes)) ;
            (FI_all_temp(3,monkey1 & cellTypes)) ;
            (FI_all_temp(3,monkey2 & cellTypes)) ;
            },...
            'CombinedIndex',15,'PlotCombinedOnly', 1, ...
            'Xlabel','FI (Vest + Vis)','Ylabel','FI (Combined)',...
            'FaceColors',{'none','none','k','k'},'Markers',{'o','^'},'EdgeColors',{'k','k','k','k'},...
            'LineStyles',{'k:','k:','k:','k:','k-'},'MarkerSize',marker_size,'SameScale',1,...
            'figN',figN,... 'XHist',20,'YHist',20,'XHistStyle','stacked','YHistStyle','stacked',
            'SameScale',1,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        
        delete([h.group(1:2).line]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');         SetFigure(20);
        axis square;
        
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot(FI_temp_vest_plus_vis(1,:),FI_all_temp(3,:),'visible','off'); hold on;
        set([gca h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, [select_tcells],1});
        
        clear Choice_pref_all_temp;
        axis tight;
        
        
    end

% 1.2. FI VS CPref ====
    function f6p5p1p2(debug)
        if debug;  dbstack;  keyboard;  end
        if isempty(fisherSimpleGu), f6p5p1(0); end
        
        monkeys = xls_num{1}(:,header.Monkey);
        monkey1 = monkeys == 15; monkey1 = monkey1(findForFisher)';
        monkey2 = monkeys == 13; monkey2 = monkey2(findForFisher)';
        
        % Which info?
        FI_all_temp = squeeze(mean(fisherSimpleGu(:,FIScatterTimeRange,:,FIScatterHeadingRange,FIScatterVariance),2))';
        tt = 1;
        
        for k = 1:3  % Plot it separately
            
            set(figure(figN),'name','FI vs. Cpref','pos',[17 514 1151 449]);
            
            cpref_sig_k = Choice_pref_p_value_all(k,findForFisher,tt) < 0.01;
            
            Choice_pref_all_temp = abs(Choice_pref_all(:,findForFisher,tt));
            
            h = LinearCorrelation({
                (Choice_pref_all_temp(k, monkey1 & cpref_sig_k ));
                (Choice_pref_all_temp(k, monkey2 & cpref_sig_k ));
                (Choice_pref_all_temp(k, monkey1 & ~cpref_sig_k)) ;
                (Choice_pref_all_temp(k, monkey2 & ~cpref_sig_k )) ;
                },...
                {
                (FI_all_temp(k,monkey1 & cpref_sig_k)) ;
                (FI_all_temp(k,monkey2 & cpref_sig_k ));
                (FI_all_temp(k,monkey1 & ~cpref_sig_k)) ;
                (FI_all_temp(k,monkey2 & ~cpref_sig_k)) ;
                },...
                'CombinedIndex',15,'PlotCombinedOnly', 1, ...
                'Ylabel','FI','Xlabel','ChoicePref',...
                'FaceColors',{condition_colors(k,:),condition_colors(k,:),'none','none'},'Markers',{'o','^'},...
                'EdgeColors',{'k','k','k','k'},...
                'LineStyles',{'k:','k:','k:','k:','k-'},'MarkerSize',marker_size,...
                'figN',figN,'XHist',20,'YHist',20,...
                'XHistStyle','stacked','YHistStyle','stacked',...
                'Method','Spearman','FittingMethod',2); figN = figN + 1;
            
            % delete([h.group(1:6).line]);
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');       %   SetFigure(20);
            
            % Annotate tcells
            %             plot((FI_all_temp(k,select_tcells(select_tcells))),(FI_all_temp(3,select_tcells(select_tcells))),...
            %                 '+','markersize',tcell_cross_size,'color','k','linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            h_line = plot(Choice_pref_all_temp(k,:),FI_all_temp(k,:),'visible','off'); hold on;
            set([gca h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, [select_tcells],1});
            
            axis tight;
        end
        %             for i = 1:3
        %                 set(h.group(i).dots,'color',colors(k,:));
        %             end
        
    end

FIDora_partial_corr_all = []; FIDora_partial_corr_all_flipped = []; FIDora_partial_corr_all_Rsquare = [];
FIDora_conditional_slope2overVar = []; FIDora_conditional_choiceslope2overVar = [];

    function f6p5p2(debug)  % Dora's partial corr / multivariate linear regression slope
        if debug;  dbstack;  keyboard;  end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        j = 1;
        select_for_partial = select_tcells;
        
        FIDora_binSize = 200; % in ms
        FIDora_stepSize = 10;
        FIDora_tCenters = FIDora_binSize/2 : FIDora_stepSize : 1500 - FIDora_binSize/2 ;
        %         conditional_variance_min_reps = 10;
        % deleted by ZZ 20210107
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        find_for_partial = find(select_for_partial);
        
        % Calculate partial correlation coefficients for each cell
        
        if isempty(FIDora_partial_corr_all)
            
            FIDora_partial_corr_all = nan(sum(select_for_partial),length(FIDora_tCenters),2,3);  % [cell number, time epoch, (heading coeff, choice coeff), stim_type]
            FIDora_partial_corr_all_flipped = nan(sum(select_for_partial),length(FIDora_tCenters),2,3);  % For plotting partial corr over time
            FIDora_conditional_slope2overVar = nan(sum(select_for_partial),length(FIDora_tCenters),1,3);
            
            progressbar('cell num');
            
            for i = 1:sum(select_for_partial)  % For cells
                
                this_raw_spike_in_bin = group_result(find_for_partial(i)).mat_raw_PSTH.spike_aligned{1,j};
                this_time = group_result(find_for_partial(i)).mat_raw_PSTH.spike_aligned{2,j};
                this_stim_type = group_result(find_for_partial(i)).mat_raw_PSTH.stim_type_per_trial;
                this_heading = group_result(find_for_partial(i)).mat_raw_PSTH.heading_per_trial;
                this_choice = group_result(find_for_partial(i)).mat_raw_PSTH.choice_per_trial;
                
                for tt = 1 : length(FIDora_tCenters)
                    
                    count_win = FIDora_tCenters(tt) - FIDora_binSize/2 <= this_time & this_time <= FIDora_tCenters(tt) + FIDora_binSize/2;
                    
                    for k = 1:3
                        if isempty(find(this_stim_type==k, 1)); continue; end
                        
                        % --- 1. Partial correlation ---
                        X=[];
                        X(:,1) = sum(this_raw_spike_in_bin(this_stim_type==k,count_win),2)...
                            / FIDora_binSize * 1e3; % Average firing rate in Hz
                        X(:,2) = this_heading(this_stim_type==k);
                        X(:,3) = this_choice(this_stim_type==k);
                        
                        r = partialcorr(X);
                        FIDora_partial_corr_all(i,tt,:,k) = r(1,2:3);
                        
                        % --- 2. Partial FI ---
                        % -     2.1 Multivariate linear regression of beta ---
                        normalizedX3 = X(:,3)./rms(X(:,3)).*rms(X(:,2)); % Normalize choice to headings' RMS (to keep sensory heading comparable with the traditional FI)
                        
                        % Orthogonalize the regressors
                        [score,~] = mgson([X(:,2) normalizedX3]);
                        
%                         [~,score] = pca([X(:,2) normalizedX3]); 
                        coeff = glmfit(score,X(:,1));
                        
%                         coeff = glmfit([X(:,2) normalizedX3],X(:,1));
                        conditional_sensory_slope = coeff(2) * (180/pi); % Turn to rad
                        conditional_choice_slope = coeff(3) * (180/pi); % Turn to rad, although it's choice
                        
                        % -     2.2 Compute conditional variance
                        unique_heading = unique(this_heading);
                        conditional_variance_matrix = nan(2,length(unique_heading));
                        
                        for hh = 1:length(unique(this_heading))
                            for cc = LEFT:RIGHT
                                this_rates = X( X(:,2) == unique_heading(hh) & X(:,3) == cc,1);
                                %                                 if length(this_rates) >= conditional_variance_min_reps
                                conditional_variance_matrix(cc,hh) = var(this_rates);
                                %                                 end
                            end
                        end
                        
                        conditional_variance_mean = nanmean(conditional_variance_matrix(:));
                        
                        FIDora_conditional_slope2overVar(i,tt,1,k) = conditional_sensory_slope^2/conditional_variance_mean;
                        FIDora_conditional_choiceslope2overVar(i,tt,1,k) = conditional_choice_slope^2/conditional_variance_mean;
                    end
                end
                
                % Use flipped partial corr (Fig.6 of Zaidel 2017)
                if group_result(find_for_partial(i)).PREF_PSTH == LEFT  % Flip according to PREF of this cell
                    FIDora_partial_corr_all_flipped(i,:,:,:) = - FIDora_partial_corr_all(i,:,:,:);
                end
                
                progressbar(i/sum(select_for_partial));
            end
            
            % Use R^2 (Fig.4 of Zaidel 2017)
            FIDora_partial_corr_all_Rsquare = FIDora_partial_corr_all.^2;
        end
        
        % ========  Plot partial correlation over time (Zaidel 2017, Figure 6)  HH20180821
        set(figure(082215),'name',sprintf('Partial correlation over time, j = %g, "any significant" out of N = %g, "%s" cells',...
            j,sum(select_for_partial),t_criterion_txt)); clf;
        set(gcf,'uni','norm','pos',[0.016       0.257       0.965       0.528]);
        
        hhcc_text = {'Heading','Choice'};
        
        %         for hhcc = 1:2
        
        for k = 1:3
            % Flip
            SeriesComparison(squeeze(FIDora_partial_corr_all_flipped(:,:,:,k)), FIDora_tCenters,...
                'SEM',1,'Errorbar',2,'Axes',subplot(2,5,k),'Transparent',transparent,'colors',[condition_colors(k,:); 0 0 0]);         legend off;
            
            % Gaussian vel
            axis tight; plot([0 0],ylim,'k--'); xlim([-100 1600]); plot(xlim,[0 0],'k--'); % ylim([0 max(ylim)*1.1])
            plot([1500 1500], ylim,'k--')
            plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2.5,'color',[0.6 0.6 0.6]);
            
            % R^2
            SeriesComparison(squeeze(FIDora_partial_corr_all_Rsquare(:,:,:,k)), FIDora_tCenters,...
                'SEM',1,'Errorbar',2,'Axes',subplot(2,5,5 + k),'Transparent',transparent,'colors',[condition_colors(k,:); 0 0 0]);         legend off;
            
            % Gaussian vel
            axis tight; plot([0 0],ylim,'k--'); xlim([-100 1600]); plot(xlim,[0 0],'k--'); % ylim([0 max(ylim)*1.1])
            plot([1500 1500], ylim,'k--')
            plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2.5,'color',[0.6 0.6 0.6]);
            
        end
        
        for hhcc = 1:2
            
            SeriesComparison(squeeze(FIDora_partial_corr_all_flipped(:,:,hhcc,:)), FIDora_tCenters,...
                'SEM',1,'Errorbar',2,'Axes',subplot(2,5,3+hhcc),'Transparent',transparent,'colors',condition_colors);     legend off;
            
            title([hhcc_text{hhcc} ' partial R flip'])
            
            % Gaussian vel
            xlim([-100 1600]); plot(xlim,[0 0],'k--'); % ylim([0 max(ylim)*1.1])
            plot([0 0],ylim,'k--'); plot([1500 1500], ylim,'k--')
            plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2.5,'color',[0.6 0.6 0.6]);
            
            
            SeriesComparison(squeeze(FIDora_partial_corr_all_Rsquare(:,:,hhcc,:)), FIDora_tCenters,...
                'SEM',1,'Errorbar',2,'Axes',subplot(2,5,5+3+hhcc),'Transparent',transparent,'colors',condition_colors);  legend off;
            title([hhcc_text{hhcc} ' partial R^2'])
            
            %         sumFisherMean = squeeze(mean(sumBoots,1));
            %         sumFisherVestPlusVIs = sumFisherMean(:,1) + sumFisherMean(:,2);
            %         plot(CP_ts{1}(plotRange),sumFisherVestPlusVIs,'m','linew',2);
            
            % Gaussian vel
            axis tight; plot([0 0],ylim,'k--'); plot([1500 1500], ylim,'k--')
            xlim([-100 1600]); plot(xlim,[0 0],'k--'); % ylim([0 max(ylim)*1.1])
            plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2.5,'color',[0.6 0.6 0.6]);
            
        end
        SetFigure(15);
        
        figure(082217); clf
        set(gcf,'uni','norm','pos',[0.053       0.422       0.605       0.387]);
        
        % -------- Plotting FI(t) ---------
        bootN = 2000;
        FIDora_conditional_slope2overVar(FIDora_conditional_slope2overVar==Inf)=nan;
        FIDora_conditional_choiceslope2overVar( FIDora_conditional_choiceslope2overVar==Inf)=nan;
        % added by ZZ 20200106
        
        % Get sum of Fisher and std by bootstrap (Gu 2010)
        h = subplot(1,2,1);
        sumBoots = bootstrp(bootN,@(x)nansum(x,1),squeeze(FIDora_conditional_slope2overVar));
        %         sumBoots = bootstrp(bootN,@(x)sum(x,1),squeeze(FIDora_conditional_slope2overVar));
        sumBoots = reshape(sumBoots,bootN,[],3);
        
        SeriesComparison(sumBoots, FIDora_tCenters,...
            'SEM',0,'Errorbar',2,'axes',h,  'Transparent',transparent,'colors',condition_colors);     legend off;
        
        sumFisherMean = squeeze(nanmean(sumBoots,1));
        %         sumFisherVestPlusVIs = sumFisherMean(:,1) + sumFisherMean(:,2);
        %         plot(FIDora_tCenters,sumFisherVestPlusVIs,'m','linew',2);
        
        % Gaussian vel
        xlim([-100 1600]); % ylim([0 5000]);
        ylim([0 max(ylim)*1.1]);
        plot([0 0],ylim,'k--'); plot([1500 1500], ylim,'k--')
        plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2,'color',[0.6 0.6 0.6]);
        legend off;
        title('Conditional FI,Sensory')
        
        % ----- Partial Choice FI ----
        
        h = subplot(1,2,2);
        sumBoots = bootstrp(bootN,@(x)nansum(x,1),squeeze(FIDora_conditional_choiceslope2overVar));
        sumBoots = reshape(sumBoots,bootN,[],3);
        
        SeriesComparison(sumBoots, FIDora_tCenters,...
            'SEM',0,'Errorbar',2,'axes',h, 'Transparent',transparent,'colors',condition_colors);     legend off;
        
        sumFisherMean = squeeze(nanmean(sumBoots,1));
        %         sumFisherVestPlusVIs = sumFisherMean(:,1) + sumFisherMean(:,2);
        %         plot(FIDora_tCenters,sumFisherVestPlusVIs,'m','linew',2);
        
        % Gaussian vel
        xlim([-100 1600]); ylim([0 max(ylim)*1.1])
        plot([0 0],ylim,'k--'); plot([1500 1500], ylim,'k--')
        plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2,'color',[0.6 0.6 0.6]);
        legend off;
        title('Conditional FI, Choice')
        SetFigure(15)
    end



% 2. Decoder of heading (not choice/modality) ====  HH20170810
    function f6p5p9(debug)
        if debug;  dbstack;  keyboard;  end
        
        % Override
        %         min_reps = 15; % Anne used 20
        %         training_reps = 25; % I artificially added some trials by bootstrapping if reps < training_reps.
        %                             % This can include all cells (as Gu suggested) while not change the performance too much. @HH20150521
        
        SVM_training_epoch = [0 1700];    % This make it comparable with PCA results
        
        find_for_SVM = find(select_for_SVM);
        
        % Regroup pseudo-trials pool
        pseudo_trial_pool_decoder = cell(length(find_for_SVM),4);
        pseudo_trial_pool_testing = cell(length(find_for_SVM),4,length(SVM_testing_tcenters));
        
        for ii = 1:length(find_for_SVM)
            
            this_cell_id = find_for_SVM(ii);
            
            % ====== Get unfiltered binary spike train data with all (correct and wrong) trials to decode heading!! ======
            % Get data
            spike_aligned_in_bin = group_result(this_cell_id).mat_raw_PSTH.spike_aligned{1,j};
            t_in_bin = group_result(this_cell_id).mat_raw_PSTH.spike_aligned{2,j};
            SVM_training_epoch_ind = min(SVM_training_epoch) <= t_in_bin & t_in_bin <= max(SVM_training_epoch);
            
            raw_this = spike_aligned_in_bin(:,SVM_training_epoch_ind);
            
            % ==== Get trial conditions ====
            stim_type_per_trial = group_result(this_cell_id).mat_raw_PSTH.stim_type_per_trial';
            heading_per_trial = group_result(this_cell_id).mat_raw_PSTH.heading_per_trial';
            unique_heading = unique(heading_per_trial);
            % choice_per_trial = group_result(this_cell_id).mat_raw_PSTH.choice_per_trial; % Not important
            % correct_or_zero_per_trial = (choice_per_trial == ((heading_per_trial>0) + 1))  |(heading_per_trial == 0); % Correct or 0 heading; % Not important
            % pref_null = [group_result(this_cell_id).PREF_PSTH LEFT+RIGHT-group_result(this_cell_id).PREF_PSTH]; % Pref goes first % Not important
            
            
            % Training pool (one certain time window)
            means_this = cellfun(@(x)mean(x(:,SVM_training_epoch_ind),2),raw_this,'uniformOutput',false);
            [pseudo_trial_pool_decoder{ii,:}] = means_this{:};
            
            % Testing pool (shifting time windows)
            for ttt = 1:length(SVM_testing_tcenters)
                test_epoch_ind = SVM_testing_tcenters(ttt) - SVM_testing_span/2 <= rate_ts{j_for_SVM} & rate_ts{j_for_SVM} <= SVM_testing_tcenters(ttt) + SVM_testing_span/2;
                means_this = cellfun(@(x)mean(x(:,test_epoch_ind),2),raw_this,'uniformOutput',false);
                
                [pool_testing{ii,:,ttt}] = means_this{:};
            end
            
        end
        
        % Find cells who cross the minimal repetitions for each condition (Anne)
        lengths = cellfun(@(x)length(x),pseudo_trial_pool_decoder);
        select_for_SVM_actual = all(lengths >= min_reps_training,2);
        pseudo_trial_pool_decoder = pseudo_trial_pool_decoder(select_for_SVM_actual,:);
        pseudo_trial_pool_testing = pseudo_trial_pool_testing(select_for_SVM_actual,:,:);
        
        lengths = lengths(select_for_SVM_actual,:);
        
        % I artificially added some trials by bootstrapping if reps < training_reps.
        % This can include all cells (as Gu suggested) while not change the performance too much. @HH20150521
        
        who_needs_add_some_trials = find(lengths < training_reps);
        
        for aa = 1:length(who_needs_add_some_trials)
            this_who = who_needs_add_some_trials(aa);
            
            % Select trials
            real_reps = length(pseudo_trial_pool_decoder{this_who});
            n_to_add = training_reps - real_reps;
            trials_to_add = randperm(real_reps,n_to_add);
            
            % Add trials
            pseudo_trial_pool_decoder{this_who} = [pseudo_trial_pool_decoder{this_who}; pseudo_trial_pool_decoder{this_who}(trials_to_add)];
            [sub1,sub2] = ind2sub(size(pseudo_trial_pool_decoder),this_who);
            
            for ttt = 1:length(SVM_testing_tcenters)
                pseudo_trial_pool_testing{sub1,sub2,ttt} = [pseudo_trial_pool_testing{sub1,sub2,ttt}; pseudo_trial_pool_testing{sub1,sub2,ttt}(trials_to_add)];
            end
        end
        
        % Recalculate lengths
        lengths = cellfun(@(x)length(x),pseudo_trial_pool_decoder);
        reps_actual_training = min(lengths);
        
        % -------- Teacher signals --------
        % For choices, 1 = Contralateral, 0 = Ipsilateral
        answers_choice = [ones(1,reps_actual_training(1)) zeros(1,reps_actual_training(2)) ones(1,reps_actual_training(3)) zeros(1,reps_actual_training(4))]';
        % For modality, 0 = Vestibular, 1 = Visual
        answers_modality = [zeros(1,reps_actual_training(1)) zeros(1,reps_actual_training(2)) ones(1,reps_actual_training(3)) ones(1,reps_actual_training(4))]';
        
        %% =============  SVM Bootstrapping ============
        % Parallel computing
        if ver_num < 2014
            if matlabpool('size') == 0 ;  matlabpool;  end
        else
            if isempty(gcp('nocreate')); parpool; end
        end
        
        training_set = cell(n_bootstrap_weights);
        tic; disp('Starting SVM training...');
        
        parfor_progress(n_bootstrap_weights);
        
        parfor nn = 1:n_bootstrap_weights
            % Randomly select trials without replacement from the pseudo-trials pool
            pseudo_trial_pool_perm = cellfun(@(x)x(randperm(size(x,1)),1),pseudo_trial_pool_decoder,'uniform',false);
            
            for conds = 1:4
                training_set{nn}{conds} = cell2mat(cellfun(@(x)x(1:reps_actual_training(conds),:),pseudo_trial_pool_perm(:,conds),'uniform',false)');
            end
            
            firing_for_this_training = cell2mat(training_set{nn}');
            
            % Train two SVM classifiers (choice and modality)
            svm_training_choice(nn) = svmtrain(firing_for_this_training,answers_choice,'box',1e-5,'tol',1e-7,'autoscale',0);
            svm_training_modality(nn) = svmtrain(firing_for_this_training,answers_modality,'box',1e-5,'tol',1e-7,'autoscale',0);
            
            parfor_progress;
        end
        parfor_progress(0);
        toc; disp('Done!');
        
        %% Averaged weights (bagging)
        %         weights_svm_choice_allbootstrap = cell2mat(arrayfun(@(x)-x.SupportVectors'*x.Alpha,svm_training_choice,'uniform',0));
        %         weights_svm_modality_allbootstrap = cell2mat(arrayfun(@(x)-x.SupportVectors'*x.Alpha,svm_training_modality,'uniform',0));
        
        % By default, 'autoscale' = true, so here I should rescale back the weight! HH20170613
        weights_svm_choice_allbootstrap = cell2mat(arrayfun(@(x)-x.SupportVectors'*x.Alpha,svm_training_choice,'uniform',0));
        weights_svm_modality_allbootstrap = cell2mat(arrayfun(@(x)-x.SupportVectors'*x.Alpha,svm_training_modality,'uniform',0));
        
        % Orthogonal test
        angles = acos(sum(weights_svm_choice_allbootstrap.*weights_svm_modality_allbootstrap,1)./...
            sqrt(sum(weights_svm_choice_allbootstrap.^2,1)./ sum(weights_svm_modality_allbootstrap.^2,1)))/pi*180;
        
        weights_svm_choice_mean = mean(weights_svm_choice_allbootstrap,2);
        weights_svm_choice_mean = weights_svm_choice_mean / max(abs(weights_svm_choice_mean));
        % svm_weights_choice_sem = std(svm_weights_choice,[],2)/sqrt(n_bootstrap);
        
        weights_svm_modality_mean = mean(weights_svm_modality_allbootstrap,2);
        weights_svm_modality_mean = weights_svm_modality_mean / max(abs(weights_svm_modality_mean));
        % svm_weights_modality_sem = std(svm_weights_modality,[],2)/sqrt(n_bootstrap);
        
        %% Find proper threshold for SVM
        training_set = cell(n_bootstrap_threshold);
        bestZ_choice = nan(n_bootstrap_threshold,1);
        bestZ_modality = nan(n_bootstrap_threshold,1);
        tic; disp('Find proper threshold for SVM decoders...');
        
        parfor nn = 1:n_bootstrap_threshold
            % Randomly select trials without replacement from the pseudo-trials pool
            pseudo_trial_pool_perm = cellfun(@(x)x(randperm(size(x,1)),1),pseudo_trial_pool_decoder,'uniform',false);
            
            for conds = 1:4
                training_set{nn}{conds} = cell2mat(cellfun(@(x)x(1:reps_actual_training(conds),:),pseudo_trial_pool_perm(:,conds),'uniform',false)');
            end
            
            proj_choice_on_choice = cell2mat(training_set{nn}') * weights_svm_choice_mean;
            proj_modality_on_modality = cell2mat(training_set{nn}') * weights_svm_modality_mean;
            
            % Find the best threshold (the most northwestern point on ROC curve)
            [~,bestZ_choice(nn)] = rocN(proj_choice_on_choice(logical(answers_choice)),proj_choice_on_choice(~logical(answers_choice)));
            [~,bestZ_modality(nn)] = rocN(proj_modality_on_modality(logical(answers_modality)),proj_modality_on_modality(~logical(answers_modality)));
            
        end
        toc; disp('Done!');
        
        thres_choice = mean(bestZ_choice);
        thres_modality = mean(bestZ_modality);
        
        % ================ SVM Training End ===============
        
    end


%%%%%%%%%%%%%%%%%%%%%%%%%  8. TDR   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weights_TDR_PCA_SVM_mean = [];
weights_TDR_PCA_SVM_allbootstrap = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function f8p1(debug)      % TDR: PCA + SVM
        if debug;  dbstack;  keyboard;  end
        
        % --------  Targeted dimensionality reduction ----------
        if isempty(weights_PCA_B_PC)
            f5p0(0);
        end
        
        if isempty(weights_svm_choice_mean)
            f6p0(0);
        end
        
        % ------ Mean weight --------
        % Denoised matrix from PCA
        D = weights_PCA_B_PC(:,1:denoised_dim) * weights_PCA_B_PC(:,1:denoised_dim)';
        
        % Project SVM weights in subspace spanned by the first denoised_dim principal components
        beta = [weights_svm_choice_mean weights_svm_modality_mean];
        beta_pca = D * beta;
        
        % Orthogonalization by QR-decomposition
        [Q,~] = qr(beta_pca);
        weights_TDR_PCA_SVM_mean = Q(:,1:2);   weights_TDR_PCA_SVM_mean(:,1) = - weights_TDR_PCA_SVM_mean(:,1);
        
        % Contribution of naive PCs
        pc_contribution =  weights_PCA_B_PC(:,1:denoised_dim)' * weights_TDR_PCA_SVM_mean;
        figure(64); bar(pc_contribution,1) ; SetFigure(); set(gca,'xtick',1:denoised_dim);
        
        % ------ Bootstrap weights --------
        weights_TDR_PCA_SVM_allbootstrap = nan(size(weights_svm_choice_allbootstrap,1),2,size(weights_svm_choice_allbootstrap,2));
        
        for bb = 1:size(weights_svm_choice_allbootstrap,2)
            
            % Project SVM weights in subspace spanned by the first denoised_dim principal components
            beta = [weights_svm_choice_allbootstrap(:,bb) weights_svm_modality_allbootstrap(:,bb)];
            beta_pca = D * beta;
            
            % Orthogonalization by QR-decomposition
            [Q,~] = qr(beta_pca);
            weights_TDR_PCA_SVM_allbootstrap(:,:,bb) = Q(:,1:2);
            
            % Normalized by sum, not 2-norm
            weights_TDR_PCA_SVM_allbootstrap(:,1,bb) = (weights_TDR_PCA_SVM_allbootstrap(:,1,bb))/sum((weights_TDR_PCA_SVM_allbootstrap(:,1,bb)),1);
            weights_TDR_PCA_SVM_allbootstrap(:,2,bb) = (weights_TDR_PCA_SVM_allbootstrap(:,2,bb))/sum((weights_TDR_PCA_SVM_allbootstrap(:,2,bb)),1);
            
        end
        
        % -------- Weights vs cell properties ----------
        Weights_Property_Correlation(weights_TDR_PCA_SVM_mean,...
            {'Weights for TDR axis 1 (choice)','Weights for TDR axis 2 (modality)'}, select_for_PCA_B);
        
    end
    function f8p2(debug)      % TDR: PCA + SVM, correct trials / angles
        if debug;  dbstack;  keyboard;  end
        
        if isempty(weights_TDR_PCA_SVM_mean)
            f8p1(0);
        end
        
        % ------ Projection of raw data on task-related axes -------
        %         h = Weighted_sum_PSTH({weights_TDR_PCA_SVM_mean(:,1),weights_TDR_PCA_SVM_mean(:,2)},{'Weighted by TDR\_PCA\_SVM\_choice','Weighted by TDR\_PCA\_SVM\_modality'},...
        %             {select_for_SVM select_for_SVM});
        
        h = Weighted_sum_PSTH({squeeze(weights_TDR_PCA_SVM_allbootstrap(:,1,:)),squeeze(weights_TDR_PCA_SVM_allbootstrap(:,2,:))},...
            {'Weighted by TDR\_PCA\_SVM\_choice','Weighted by TDR\_PCA\_SVM\_modality'},...
            {select_for_SVM select_for_SVM});
        
        %
        figure(904);  clf; hold on;
        for k = 1:3
            %             plot(rate_ts{1}(1:end-1) ,diff(h.projected_all_mean{2}(k,:)-h.projected_all_mean{1}(k,:)),'color',colors(k,:),'linew',2);
            %             plot(rate_ts{1}(1:end-1) ,diff(h.projected_all_mean{1}(k*2-1,:)-h.projected_all_mean{1}(k*2,:)),'color',colors(k,:),'linew',2);
            plot(rate_ts{1}(1:end) ,h.projected_all_mean{1}{1}(k,:),'color',condition_colors(k,:),'linew',2);
        end
        SetFigure;
        
        %% ------ Tuning curves ------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        j = 1;
        tuning_centers = [ mean(time_markers{j}(1,1:2)) + 150  % Stim center + sensory delay
            mean(time_markers{j}(1,3)) - group_result(representative_cell).mat_raw_PSTH.binSize_CP/2   % Pre-sac epoch
            mean(time_markers{j}(1,3)) + group_result(representative_cell).mat_raw_PSTH.binSize_CP/2];  % Post-sac epoch
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [~,center_t_ind] = min(abs(rate_ts{j} - tuning_centers(1)));
        [~,pre_t_ind] = min(abs(rate_ts{j} - tuning_centers(2)));
        [~,post_t_ind] = min(abs(rate_ts{j} - tuning_centers(3)));
        
        tuning_time_phase = [center_t_ind pre_t_ind post_t_ind];
        tuning_time_phase_title = {'Stimulus', 'Pre-sac', 'Post-sac'};
        
        unique_heading = group_result(representative_cell).mat_raw_PSTH.CP{j,1}.raw_CP_result{1}.Neu_tuning(:,1);
        unique_heading = [unique_heading(unique_heading<0);-0.00001; 0.00001; unique_heading(unique_heading>0)]; % Add another 0 heading
        
        % Get means and sems
        %         tuning_mean_correctonly = permute(h.projected_angles_mean{j},[2 1 3]);
        tuning_mean_correctonly = h.projected_angles_mean{j};
        tuning_sem_correctonly = h.projected_angles_sem{j};
        
        smooth_time_win = 500;
        for hh = 1:size(tuning_mean_correctonly,1)
            for kk = 1:size(tuning_mean_correctonly,3)
                tuning_mean_correctonly(hh,:,kk) = smooth(tuning_mean_correctonly(hh,:,kk), round(smooth_time_win/mean(diff(rate_ts{1}))));
            end
        end
        
        % Sort to match 'unique_heading'
        tuning_mean_correctonly =  [tuning_mean_correctonly(end:-2:1,:,:); tuning_mean_correctonly(1:2:end,:,:)];
        
        % Plotting tuning curves at three time points
        
        set(figure(2099+figN),'pos',[9 546 1378 414],'name',['Tuning curve, j = ' num2str(j)]); clf; figN = figN+1;
        
        for pp = 1:length(tuning_time_phase)
            
            for k = 1:3  % For each stim type
                
                % Plotting
                subplot(1,length(tuning_time_phase),pp );  hold on; ylabel('Correct only');
                plot(unique_heading(end/2+1:end),tuning_mean_correctonly(end/2+1:end,tuning_time_phase(pp),k),'-o','markersize',9,'color',condition_colors(k,:),'markerfacecol',condition_colors(k,:),'LineWid',2);
                plot(unique_heading(1:end/2),tuning_mean_correctonly(1:end/2,tuning_time_phase(pp),k),'-o','markersize',9,'color',condition_colors(k,:),'markerfacecol','none','LineWid',2);
                
                h_tmp = errorbar(unique_heading,tuning_mean_correctonly(:,tuning_time_phase(pp),k),tuning_sem_correctonly(:,tuning_time_phase(pp),k),'color',condition_colors(k,:),'LineWid',2);
                errorbar_tick(h_tmp,10000);
                
                title(num2str(tuning_centers(pp)));
                axis tight;
                
                SetFigure(20);
            end
        end
        
        
        %% Tuning: Hotgram
        
        to_evolve =  {tuning_mean_correctonly};
        to_evolve_title = {'Correct only'};
        
        for co = 1
            set(figure(2099+figN),'pos',[82 99 927 855],'name','Evolve of tuning curves'); clf; figN = figN+1;
            
            [X,Y] = meshgrid(rate_ts{j},unique_heading);
            [Xq,Yq] = meshgrid(linspace(min(rate_ts{j}),max(rate_ts{j}),1000),linspace(min(unique_heading),max(unique_heading),100));
            
            
            for k = 1:3  % For each stim type
                subplot_tight(1,3,k,0.02,0.1);
                %                 surf(Xq,Yq,interp2(X,Y,to_evolve{co}(:,:,k),Xq,Yq),'edgecolor','none'); view(90,-90); hold on;
                contourf(Xq,Yq,interp2(X,Y,to_evolve{co}(:,:,k),Xq,Yq),40,'edgecolor','k'); view(90,-90); hold on;
                
                if k>1 ;set(gca,'xtick',[]); end
                if k==2 ; title(to_evolve_title{co}); end
                
                caxis([min(to_evolve{co}(:)) 1.05*max(to_evolve{co}(:))]); % Use the same color range
                xlim([min(CP_ts{j}) max(CP_ts{j})]); ylim([min(unique_heading) max(unique_heading)]);
                
                for tt = 1:3
                    plot3([1 1] * time_markers{j}(1,tt),ylim,-3*[1 1],'k','linestyle',marker_for_time_markers{j}{tt},'linew',2);
                end
                
                load('colormap_for_tuning_hotgram_TDR','colormap_for_tuning_hotgram_TDR');
                colormap(colormap_for_tuning_hotgram_TDR);
                
                %                 axis tight;
                % Gaussian vel
                plot(Gauss_vel(:,1) + time_markers{j}(1),Gauss_vel(:,2)*range(ylim)/2 + min(ylim),'--','linew',2.5,'color',[0.6 0.6 0.6]);
                
            end
            
            SetFigure();
        end
        
    end


    function f9p1(debug)      % Cell counter
        if debug;  dbstack;  keyboard;  end
        
        % Plot cell_num v.s. experimental days to see if I've been too lazy (which is beyond doubt)...
        dates = xls_num{1}(:,header.Date);
        day_counter = datenum(num2str(dates),'yyyymmdd');
        
        set(figure(9999),'name','My Progress','position',[17 258 649 693]); clf;
        plot([day_counter; today],[cumsum(select_all); sum(select_all)],'k-','linew',3); hold on;
        plot([day_counter; today],[cumsum(select_sus); sum(select_sus)],'r-','linew',3);
        plot([day_counter; today],[cumsum(select_tcells); sum(select_tcells)],'g-','linew',3);
        legend({'SUs + MUs','SUs (bottom line)','T SUs'},'location','best');
        
        plot(xlim,[100 100],'r--');
        plot([today today],ylim,'k-');
        
        xlims = xlim;
        set(gca,'XTick',linspace(xlims(1),xlims(2),12));
        datetick('x',26,'KeepTicks');
        
        rotateXLabels(gca,30);
        axis tight;
        SetFigure();
        
        % After running this part of code, I am literally crying...
        % Come on man!!!!!!!
        
        % Plot prediction ratio
        set(figure(9998),'name','Prediction ratio','position',[686 539 900 412]); clf;
        h_line = plot(Psy_pred_ratio,'ko','markerfacecol','k','markersize',10); hold on;
        
        plot(find(select_tcells),Psy_pred_ratio(select_tcells),'+r','linew',2,'markersize',10);
        plot(xlim,[1 1],'k--'); ylabel('Prediction ratio'); set(gca,'yscale','log','ytick',[0.6 1 2]);
        SetFigure(15);
        
        % Plot threshold for single cues
        %          Psy_pred_ratio_vestibular_visual = sort(Psy_pred_ratio_vestibular_visual,2);
        plot(Psy_pred_ratio_vestibular_visual(:,1) ,'b-','linew',2);
        plot(Psy_pred_ratio_vestibular_visual(:,2) ,'r-','linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        set([gca h_line],'ButtonDownFcn',{@Show_individual_cell,h_line,true(size(Psy_pred_ratio))});
        
        % --- Venn diagram for cell counters -- HH20160915
        % Fig.A: Choice preferences for three modalities
        venn_HD_1 = select_bottom_line_all_monkey & group_ChoicePreference_pvalue(1,:,1)' < 0.01;
        venn_HD_2 = select_bottom_line_all_monkey & group_ChoicePreference_pvalue(2,:,1)' < 0.01;
        venn_HD_3 = select_bottom_line_all_monkey & group_ChoicePreference_pvalue(3,:,1)' < 0.01;
        venn_HD_12 = venn_HD_1 & venn_HD_2;
        venn_HD_13 = venn_HD_1 & venn_HD_3;
        venn_HD_23 = venn_HD_2 & venn_HD_3;
        venn_HD_123 = venn_HD_1 & venn_HD_2 & venn_HD_3;
        
        fprintf('HD_1, HD_2, HD_3: %s\n', num2str([sum(venn_HD_1),sum(venn_HD_2),sum(venn_HD_3),...
            sum(venn_HD_12),sum(venn_HD_13),sum(venn_HD_23),sum(venn_HD_123)]));
        
        figure(9997);  clf
        [H,S]=venn([sum(venn_HD_1),sum(venn_HD_2),sum(venn_HD_3)],...
            [sum(venn_HD_12),sum(venn_HD_13),sum(venn_HD_23),sum(venn_HD_123)],...
            'FaceColor',{condition_colors(1,:),condition_colors(2,:),condition_colors(3,:)},'FaceAlpha',0.5,'Edgecolor',{condition_colors(1,:),condition_colors(2,:),condition_colors(3,:)});
        for i = 1:length(S.ZoneCentroid)
            text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), num2str(S.ZonePop(i)));
        end
        title(sprintf('HD_1, HD_2, HD_3: %s, total = %g',num2str([S.CirclePop S.IntersectPop]),sum(S.ZonePop)));
        SetFigure();
        axis off;
        
        % Fig.B: Memsac and HD
        venn_HD_any = venn_HD_1 | venn_HD_2 | venn_HD_3;
        venn_HD_all = venn_HD_1 & venn_HD_2 & venn_HD_3;
        venn_Mem = select_bottom_line_all_monkey & (group_MemSac_ps(:,3) < 0.05 | group_TwoPtMemSac_ps(:,3) < 0.05);
        venn_HD_any_all = venn_HD_any & venn_HD_all;
        venn_HD_any_Mem = venn_HD_any & venn_Mem;
        venn_HD_all_Mem = venn_HD_all & venn_Mem;
        venn_HD_any_all_Mem = venn_HD_any & venn_HD_all & venn_Mem;
        
        fprintf('HD_any, HD_all, Memsac: %s\n', num2str([sum(venn_HD_any),sum(venn_HD_all),sum(venn_Mem),...
            sum(venn_HD_any_all),sum(venn_HD_any_Mem),sum(venn_HD_all_Mem),sum(venn_HD_any_all_Mem)]));
        fprintf('HD_any, Memsac: %s\n', num2str([sum(venn_HD_any),sum(venn_Mem),...
            sum(venn_HD_any_Mem)]));
        
        % Fig.C: Memsac and HD (simple)
        figure(9995);  clf
        [H,S]=venn([sum(venn_HD_any),sum(venn_Mem)],...
            [sum(venn_HD_any_Mem)],...
            'FaceColor',{'none','none'},'Edgecolor',{condition_colors(1,:),condition_colors(2,:)});
        for i = 1:length(S.ZoneCentroid)
            text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), num2str(S.ZonePop(i)));
        end
        title(sprintf('HD_any, Memsac: %s, total = %g',num2str([S.CirclePop S.IntersectPop]),sum(S.ZonePop)));
        
        
    end
    function f9p2(debug)      % Others 1. Targ first vs Target Last (All SU + MU)
        if debug;  dbstack;  keyboard;  end
        
        %%
        
        j = 1;
        selectCells = (xls_num{1}(:,header.Units_RealSU) >=0 ) & (xls_num{1}(:,header.HD_TargFirst)~=0) & (~isnan(ChoiceDiv_All{j}(:,1,k)));
        
        figure(2099+figN); clf; figN = figN+1;
        
        for k = 1:3
            ys = nanmean(ChoiceDiv_All{j}(selectCells,:,k));
            errors = nanstd(ChoiceDiv_All{j}(selectCells,:,k))/sqrt(sum(selectCells));
            %             h1 = shadedErrorBar(rate_ts{j},ys,errors * 0,{'Color',colors(k,:)},transparent);
            
            h1.mainLine  = plot(rate_ts{j},ys,'Color',condition_colors(k,:));
            set(h1.mainLine,'LineWidth',3);
            
            hold on;
        end
        
        str1 = sprintf('Choice divergence (All MU & SU, N = %g)',sum(selectCells));
        
        %%% Choice divergence: Targ last
        
        selectCells = (xls_num{1}(:,header.HD_TargFirst)==0) & (~isnan(ChoiceDiv_All{j}(:,1,k)));
        
        % figure(2099+figN); clf; figN = figN+1;
        
        for k = 1:3
            ys = mean(ChoiceDiv_All{j}(selectCells,:,k));
            errors = std(ChoiceDiv_All{j}(selectCells,:,k))/sqrt(sum(selectCells));
            
            %             h2 = shadedErrorBar(rate_ts{j},ys,errors * 0,{'Linestyle','-.','Color',colors(k,:) * 0.5 + [.5 .5 .5] },transparent);
            
            h2.mainLine = plot(rate_ts{j},ys,'Linestyle',':','Color',condition_colors(k,:) );
            set(h2.mainLine,'LineWidth',3)
            hold on;
        end
        
        str2 = sprintf('Choice divergence (Target Last, N = %g)',sum(selectCells));
        axis tight;
        
        for tt = 1:3
            plot([1 1] * time_markers{j}(1,tt),ylim,'k','linestyle',marker_for_time_markers{j}{tt},'linew',1.5);
        end
        
        SetFigure();
        
        l = legend([h1.mainLine h2.mainLine],str1,str2);
        set(l,'FontSize',10,'box','off','location','best');
        
        plot(xlim, [0 0],'k--');
        

    end

    function f9p9(debug)      % Export Associated Memsac Files
        if debug;  dbstack;  keyboard;  end
        
        % Read batch file with all cells
        f = fopen('Z:\Data\Tempo\Batch\20150725_BP_allAreas_m5_m10.m');
        line = fgetl(f);  ii = 0;
        
        while line ~= -1
            if line(1) ~= '%'
                ii = ii + 1;
                all_ori_batch{ii} = line;
                tmp = textscan(line,'%s');
                all_ori_name{ii} = tmp{1}{2};
                all_ori_spikeChan(ii) = str2double(tmp{1}{end});
            end
            line = fgetl(f);
        end
        fclose(f);
        
        % HH20180608
        exportPath = 'Z:\Data\Tempo\Batch\ExportedAssociatedMemsacFiles\';
        new_f = fopen([exportPath 'batch_file.m'],'a');
        
        progressbar;
        for ii = 1:length(group_result)
            
            PATH = 'Z:\Data\MOOG\';
            
            if xls_num{1}(ii,header.Monkey) == 5
                monkeyName = 'Polo';
            else
                monkeyName = 'Messi';
            end
            
            if ~isempty(group_result(ii).mat_raw_MemSac)
                FILE = group_result(ii).mat_raw_MemSac.FILE;
                SpikeChan = group_result(ii).mat_raw_MemSac.SpikeChan;
                %                 htbOK = copyfile([PATH monkeyName '\raw\' FILE '.htb'],[exportPath FILE '.htb']);
                %                 logOK= copyfile([PATH monkeyName '\raw\' FILE '.log'],[exportPath FILE '.log']);
                %                 matOK= copyfile([PATH monkeyName '\Analysis\SortedSpikes2\' FILE '.mat'],[exportPath FILE '.mat']);
                
                % Write new batch file
                where = find(strcmp(all_ori_name,FILE) & all_ori_spikeChan == SpikeChan);
                if length(where) ~= 1 % Should be unique
                    error('No batch info?')
                else
                    fprintf(new_f,'%s\n',all_ori_batch{where});
                end
            end
            
            if ~isempty(group_result(ii).mat_raw_2ptMemSac)
                FILE = group_result(ii).mat_raw_2ptMemSac.FILE;
                SpikeChan = group_result(ii).mat_raw_2ptMemSac.SpikeChan;
                %                 htbOK = copyfile([PATH monkeyName '\raw\' FILE '.htb'],[exportPath FILE '.htb']);
                %                 logOK= copyfile([PATH monkeyName '\raw\' FILE '.log'],[exportPath FILE '.log']);
                %                 matOK= copyfile([PATH monkeyName '\Analysis\SortedSpikes2\' FILE '.mat'],[exportPath FILE '.mat']);
                
                % Write new batch file
                where = find(strcmp(all_ori_name,FILE) & all_ori_spikeChan == SpikeChan);
                if length(where) ~= 1 % Should be unique
                    error('No batch info?')
                else
                    fprintf(new_f,'%s\n',all_ori_batch{where});
                end
            end
            
            
            progressbar(ii/length(group_result));
        end
        
        fclose(new_f);
        
        edit([exportPath 'batch_file.m']);
    end


    function Weights_Property_Correlation(weights,txt,select)  % Compare different population weights with cell properties
        % weights: [choice weight, modality weight]
        
        % ===========   1. Choice vs modality weights ============
        
        if size(weights,2) == 2
            h = LinearCorrelation({
                (weights(:,1));
                },...
                {
                weights(:,2) ;
                },...
                'CombinedIndex',[],...
                'Ylabel',txt{2},'Xlabel',txt{1},...
                'FaceColors',{'k'},'Markers',{'o'},...
                'LineStyles',{'k-'},'MarkerSize',12,...
                'figN',figN,'XHist',20,'YHist',20,...
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Pearson','FittingMethod',2); figN = figN + 1;
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
            
            % Annotate tcells
            plot(weights(select_tcells(select),1),weights(select_tcells(select),2),'+','markersize',16,'color',condition_colors(2,:),'linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            set([gca h.group.dots],'ButtonDownFcn',{@Show_individual_cell,h.group.dots,select});
            
        end
        
        % ===========   2. Correlations =============
        % ---------  Choice weights and Choice preference
        tt = 1;
        k = 1;
        
        set(figure(figN),'name','Choice weights and Choice preference','pos',[17 514 1151 449]);
        
        cpref_sig = Choice_pref_p_value_all(k,:,tt) < 0.01;
        
        h = LinearCorrelation({
            (weights(~cpref_sig,1));
            (weights(cpref_sig,1));
            },...
            {
            abs(Choice_pref_all(k,~cpref_sig ,tt)) ;
            abs(Choice_pref_all(k,cpref_sig,tt)) ;
            },...
            'CombinedIndex',[3],...
            'Ylabel','abs(Choice preference)','Xlabel',txt{1},...
            'FaceColors',{'none','k'},'Markers',{'o'},'EdgeColors',{'k'},...
            'LineStyles',{'k-'},'MarkerSize',12,...
            'figN',figN,'XHist',20,'YHist',20,...
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        delete([h.group(1:2).line]);
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        
        % Annotate tcells
        plot(weights(select_tcells(select),1),abs(Choice_pref_all(k,select_tcells(select),tt)),'+','markersize',16,'color',condition_colors(2,:),'linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot(weights(:,1),abs(Choice_pref_all(k,:,tt)),'visible','off'); hold on;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select_cpref_mpref});
        
        % ---------  Memsac DDI and Choice Weight
        monkeys = xls_num{1}(:,header.Monkey);  % Temporarily here. HH20150914
        monkey1 = monkeys == 15;
        monkey2 = monkeys == 13;
        
        h = LinearCorrelation({
            MemSac_indicator(select & monkey1);
            MemSac_indicator(select & monkey2);
            },...
            {
            (weights(monkey1(select),1));
            (weights(monkey2(select),1));
            },...
            'CombinedIndex',[1 2 3],'PlotCombinedOnly',1,...     % changed from [3] to [] by ZZ 20201230
            'Xlabel',MemSac_indicator_txt,'Ylabel',txt{1},...
            'FaceColors',{'none','k'},'Markers',{'o','^'},'EdgeColors',{'k'},...
            'LineStyles',{'k-'},'MarkerSize',12,...
            'figN',figN,'XHist',20,'YHist',20,...
            'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Pearson','FittingMethod',2); figN = figN + 1;
        %         delete([h.group(1:2).line]);
        
        plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
        
        % Annotate tcells
        plot(MemSac_indicator(select_tcells),weights(select_tcells(select),1),...
            '+','markersize',16,'color',condition_colors(2,:),'linew',2);
        
        % Show individual cell selected from the figure. HH20150424
        h_line = plot(MemSac_indicator(select),(weights(:,1)),'visible','off'); hold on;
        set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select});
        
        % ------------ Modality preference and Modality weight
        
        if size(weights,2) == 2
            
            tt = 1;
            k = 1; % Vis-vest
            
            set(figure(figN),'name','Modality weight vs Modality preference','pos',[17 514 1151 449]);
            
            mpref_sig = Modality_pref_p_value_all(1,:,tt)' < 0.05;
            
            h = LinearCorrelation({
                (Modality_pref_all(k,~mpref_sig & monkey1(select),tt)) ;
                (Modality_pref_all(k,~mpref_sig & monkey2(select),tt)) ;
                (Modality_pref_all(k,mpref_sig & monkey1(select),tt)) ;
                (Modality_pref_all(k,mpref_sig & monkey2(select),tt)) ;
                },...
                {
                (weights(~mpref_sig & monkey1(select),2));
                (weights(~mpref_sig & monkey2(select),2));
                (weights(mpref_sig & monkey1(select),2));
                (weights(mpref_sig & monkey2(select),2));
                },...
                'CombinedIndex',[5 10 15],'PlotCombinedOnly',1,...
                'Xlabel','Modality preference (vis - vest)','Ylabel',txt{2},...
                'FaceColors',{'none','none','k','k'},'Markers',{'o','^'},'EdgeColors',{'k'},...
                'LineStyles',{'k-'},'MarkerSize',12,...
                'figN',figN,'XHist',20,'YHist',20,...
                'XHistStyle','stacked','YHistStyle','stacked','SameScale',0,'Method','Pearson','FittingMethod',2); figN = figN + 1;
            
            %             delete([h.group(1:4).line]);
            plot(xlim,[0 0],'k--'); plot([0 0],ylim,'k--');    SetFigure(20);
            
            % Annotate tcells
            plot((Modality_pref_all(k,select_tcells(select),tt)), weights(select_tcells(select),2),...
                '+','markersize',16,'color',condition_colors(2,:),'linew',2);
            
            % Show individual cell selected from the figure. HH20150424
            h_line = plot((Modality_pref_all(k,:,tt)),weights(:,2),'visible','off'); hold on;
            set([gca [h.group.dots] h_line],'ButtonDownFcn',{@Show_individual_cell, h_line, select});
        end
        
    end


    function h = Weighted_sum_PSTH(weights,txt,selects)     % Plot weighted sum of firing rates
        
        % ------  PSTH correct only, all choices -----------
        
        if ~iscell(weights)  ;  weights = {weights};     end
        if ~iscell(selects) ;  selects = {selects};   end
        
        set(figure(figN),'pos',[34   83  584*length(weights)  504*2]); clf; figN = figN + 1;
        
        for pp = 1:length(weights)
            
            % Now this function can receive all bootstrap weights for computing error bars. @HH20150525
            nbootstrap = size(weights{1},2);
            
            %             for j = 1:2
            for j = 1:3
                
                PSTH_projected{j} = nan(nbootstrap,size(PSTH_all_raw{j},2),size(PSTH_all_raw{j},3));
                
                for bb = 1:nbootstrap
                    proj_this = weights{pp}(:,bb) / sum(weights{pp}(:,bb)); % Normalization to let the firing rate make sense
                    
                    % Projections
                    for kk = 1:size(PSTH_all_raw{j},3)
                        PSTH_projected{j}(bb,:,kk) = proj_this'*(PSTH_all_raw{j}(logical(selects{pp}),:,kk));
                    end
                end
            end
            
            
            % Note here the errorbars should be STD instead of SEM. (Not independent sampling, but bootstrapping)
            h_subplot = subplot(2,length(weights),pp);
            h_series = SeriesComparison({PSTH_projected{1} PSTH_projected{2} PSTH_projected{3}},...
                {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
                'Colors',{condition_colors(1,:),condition_colors(1,:),condition_colors(2,:),condition_colors(2,:),condition_colors(3,:),condition_colors(3,:)},'LineStyles',{'-','--','-','--','-','--'},...
                'SEM',0,'ErrorBar',0,'Xlabel',[],'Ylabel','Weighted sum of firing','axes',h_subplot);
            hold on;    legend off;     axis tight;
            xlim([min(xlim) + 200 max(xlim)-300]);   % ylim([0.1 .7]);
            
            h_subplot = subplot(2,length(weights),pp+length(weights));
            h_series = SeriesComparison({PSTH_projected{1}(:,:,1:2:end)-PSTH_projected{1}(:,:,2:2:end) PSTH_projected{2}(:,:,1:2:end)-PSTH_projected{2}(:,:,2:2:end) PSTH_projected{3}(:,:,1:2:end)-PSTH_projected{3}(:,:,2:2:end)},...
                {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
                'Colors',{condition_colors(1,:),condition_colors(2,:),condition_colors(3,:)},'LineStyles',{'-'},...
                'SEM',0,'ErrorBar',0,'Xlabel',[],'Ylabel','Weighted sum of firing','axes',h_subplot);
            hold on;    legend off;     axis tight;
            
            xlim([min(xlim) + 200 max(xlim)-300]);   % ylim([0.1 .7]);
            
            
            % Gaussian vel
            plot(Gauss_vel(:,1) + time_markers{1}(1), min(ylim) + Gauss_vel(:,2) * range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            
            legend off;
            h.fig = gcf;
            h.projected_all_allbootstrap{pp} = PSTH_projected;
            h.projected_all_mean{pp} = h_series.means;
            
            title([txt{1} ', n = ' num2str(sum(selects{1}))]);
            
        end
        
        SetFigure(15);
        
        % ---------- Different angles ----------
        
        colors_angles = colormap(gray);
        colors_angles = colors_angles(round(linspace(length(colors_angles)-15,1,5)),:);
        colors_angles = reshape(repmat(colors_angles,1,2)',3,[])';
        colors_angles = mat2cell(colors_angles,ones(10,1));
        
        set(figure(figN),'name',[txt{1} ', n = ' num2str(sum(selects{1}))],'pos',[26 78 888 879]); clf; figN = figN + 1;
        h_subplot = tight_subplot(3,2,[0.05 0.1],[0.05 0.15],[0.12 0.03]);
        
        for k = 1:3
            for j = 1:3
                
                %             for j = 1:2
                
                
                h.projected_angles_allbootstrap{j} = nan(nbootstrap,size(PSTH_correct_angles_raw{j},2),size(PSTH_correct_angles_raw{j},3),size(PSTH_correct_angles_raw{j},4));
                h.projected_angles_diff{j} = nan(nbootstrap,size(PSTH_correct_angles_raw{j},2),size(PSTH_correct_angles_raw{j},3)/2,size(PSTH_correct_angles_raw{j},4));
                
                for bb = 1:nbootstrap
                    
                    % ------- Different angles -------
                    % Deal with nans of 0 heading for some cells
                    PSTH_correct_angles_raw_this = PSTH_correct_angles_raw{j}(selects{1},:,:,k);
                    yyy = nan(size(PSTH_correct_angles_raw_this,2),size(PSTH_correct_angles_raw_this,3));
                    
                    % HH20170808
                    for hh = 1:size(yyy,2) % For each heading
                        non_nan_this_heading = ~(any(isnan(PSTH_correct_angles_raw_this(:,:,hh)),2));
                        yyy(:,hh) = PSTH_correct_angles_raw_this(non_nan_this_heading,:,hh)' ...
                            * weights{1}(non_nan_this_heading,bb) / sum(weights{1}(non_nan_this_heading,bb));
                    end
                    

                    
                    h.projected_angles_allbootstrap{j}(bb,:,:,k) = yyy;
                    
                    % -------- Difference ---------
                    yyy_diff = yyy(:,1:2:end) - yyy(:,2:2:end);
                    h.projected_angles_diff{j}(bb,:,:,k) = yyy_diff;
                    
                end
            end
            
            
            h_series = SeriesComparison({h.projected_angles_allbootstrap{1}(:,:,:,k) h.projected_angles_allbootstrap{2}(:,:,:,k) h.projected_angles_allbootstrap{3}(:,:,:,k)},...
                {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
                'Colors',colors_angles,'LineStyles',{'-','--'},...
                'SEM',0,'ErrorBar',2,'Xlabel',[],'Ylabel','Norm firing','axes',h_subplot(k));
            
            h.projected_angles_mean{1}(:,:,k) = h_series.means{1};
            h.projected_angles_mean{2}(:,:,k) = h_series.means{2};
            h.projected_angles_sem{1}(:,:,k) = h_series.errors{1};
            h.projected_angles_sem{2}(:,:,k) = h_series.errors{2};
            
            if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
            if k == 1
                title([txt{1} ', n = ' num2str(sum(selects{1}))]);
            end
            
            % Gaussian vel
            axis tight;
            plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            xlim([min(xlim) + 200 max(xlim)-300]);   % ylim([0.1 .7]);
            legend off;
            
            % ------ Plotting diff ------
            
            SeriesComparison({h.projected_angles_diff{1}(:,:,:,k) h.projected_angles_diff{2}(:,:,:,k) h.projected_angles_diff{3}(:,:,:,k)},...
                {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
                'Colors',colors_angles(1:2:end),'LineStyles',{'-'},...
                'SEM',0,'ErrorBar',2,'Xlabel',[],'Ylabel','Diff','axes',h_subplot(k+3));
            
            if k < 3 ;set(gca,'xtick',[]); else xlabel('Time (ms)'); end
            
            legend off;
            
            axis tight;
            
            % Gaussian vel
            axis tight;
            plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
            xlim([min(xlim) + 200 max(xlim)-300]);   % ylim([0.1 .7]);
            
        end
        SetFigure(15);
        
    end
    function Show_individual_cell(~,~,h_line, select_for_this, couples)    % Show individual cell selected from the figure. @HH20150424
        
        if nargin < 5
            couples = 1; % HH20150515 Show coupled cells in the figure
        end
        
        h_marker = guidata(gcbo);
        
        allX = get(h_line,'xData');
        allY = get(h_line,'yData');
        n_single_group = length(allX)/couples;
        
        % ------- Recover the cell number -------
        
        if ismember('control',get(gcf,'currentModifier'))  % Select cell from file name and show in the figure. @HH20150527
            fileN = input('Which cell do you want from the figure?    ','s');
            available_cells = find(select_for_this);
            
            if fileN(1) == '#'  % HH20160905. Direct input the original cell #
                ori_cell_no = str2double(fileN(2:end));
                ind = sum(select_for_this(1:ori_cell_no));
            else
                
                nn = 1; iffind = [];
                while nn <= length(available_cells) && isempty(iffind) % Find the first match
                    iffind = strfind(group_result(available_cells(nn)).cellID{1}{1},fileN);
                    nn = nn + 1;
                end
                
                if ~isempty(iffind)
                    ind = nn - 1;
                else
                    fprintf('Are you kidding...?\n');
                    return
                end
                
                ori_cell_no = find(cumsum(select_for_this)==ind,1);
            end
        else   % Select cell from figure
            pos = get(gca,'currentPoint'); posX = pos(1,1); posY = pos(1,2);
            [min_dis,ind] = min(abs(((posX-allX)/range(xlim)).^2+((posY-allY)/range(ylim)).^2));
            if min_dis > (range(xlim)^2+range(ylim)^2)/100 +inf ; return; end
            
            ind = mod(ind-1,n_single_group) + 1; % Deal with coupled cells
            ori_cell_no = find(cumsum(select_for_this)==ind,1);
        end
        
        % Plotting
        if ~isempty(h_marker) ; try delete(h_marker); catch ; end ;end
        
        all_inds = mod(1:length(allX),n_single_group) == mod(ind,n_single_group); % Deal with coupled cells
        h_marker = plot(allX(all_inds),allY(all_inds),'x','color','m','markersize',15,'linew',3);
        
        % Do plot
        Plot_HD([],[],ori_cell_no);
        
        guidata(gcbo,h_marker);
    end

    function Plot_HD(~,~,ori_cell_no,h_1463_axes)    % Plot associated HD traces @HH20150426
        
        
        % ------ PSTH in HD ------
        j_this = 1;
        
        for j = 1:3
            %         for j = 1:2
            ys_this{j} = group_result(ori_cell_no).mat_raw_PSTH.PSTH{j,1,1}.ys';
            sem_this{j} = group_result(ori_cell_no).mat_raw_PSTH.PSTH{j,1,1}.sem';
            ps_this{j} = group_result(ori_cell_no).mat_raw_PSTH.PSTH{j,1,1}.ps';
        end
        
        % ------ Adapted for batch plot example cells ---
        if nargin < 4
            figure(1463); clf ;
            set(gcf,'uni','norm','pos',[0.005       0.597       0.992       0.304]);
            h_1463_axes = tight_subplot(1,5,[0.1 0.06],[0.2 0.15],[0.05 0.05]);
        end
        
        % --- 1. Raw PSTHs ---
        % h_1463_PSTH = subplot_tight(1,5,1,[0.1 0.06],[0.2 0.15 0.05 0.05]);
        h_1463_PSTH = h_1463_axes(1);
        axes(h_1463_PSTH)
        
        SeriesComparison({shiftdim(ys_this{1},-1) shiftdim(ys_this{2},-1) shiftdim(ys_this{3},-1)},...
            {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
            'OverrideError',{sem_this{1}, sem_this{2} sem_this{3}},...
            'OverridePs',{ps_this{1}, ps_this{2},ps_this{3}},'ErrorBar',6,...
            'CompareIndex',[1 3 5;2 4 6],'CompareColor',{condition_colors(1,:),condition_colors(2,:),condition_colors(3,:)},...
            'Colors',{condition_colors(1,:),condition_colors(1,:),condition_colors(2,:),condition_colors(2,:),condition_colors(3,:),condition_colors(3,:)},...
            'Transparent',transparent,'LineStyles',{'-','--'},'axes',h_1463_PSTH, 'PCritical', 0.01);
        
        axis tight; legend off;
        %         xlim([-300 2300]);
        %         xticks([0 1000 2050]); xticklabels({'Stim','1000 ms','Sac'});
        
        
        id_this = group_result(ori_cell_no).cellID{1}{1};
        name_this = id_this(strfind(id_this,'u'):end);
        pos_this = group_result(ori_cell_no).position;
        
        hemi = 'L' * (pos_this(2) == 1) + 'R' * (pos_this(2) == 2);
        
        %         title(sprintf('#%g, %s, %s\n%s: AP%gML%.2gD%g',...
        %             ori_cell_no,celltype, name_this, hemi, pos_this(3:5)));
        title(sprintf('#%g, %s, %s\n%s: AP%gML%.2gD%g',...
            ori_cell_no, name_this, hemi, pos_this(3:5)));
        % changed by ZZ 20201229
        
        % xlabel('Time to saccade onset (ms)');      ylabel('Firing rate');
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{1}(1),0 + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
        
        
        % --- 2. PSTH_difference ---
        h_1463_PSTH_diff = h_1463_axes(2);
        axes(h_1463_PSTH_diff)
        
        %         for j = 1:2
        for j = 1:3
            
            psth_diff{j} = ys_this{j}(:,[1 3 5]) - ys_this{j}(:,[2 4 6]);
            psth_diff_sem{j} = sqrt(sem_this{j}(:,[1 3 5]).^2 + sem_this{j}(:,[2 4 6]));
        end
        
        SeriesComparison({shiftdim(psth_diff{1},-1) shiftdim(psth_diff{2},-1) shiftdim(psth_diff{3},-1)},...
            {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
            'OverrideError',{psth_diff_sem{1}, psth_diff_sem{2} psth_diff_sem{3}},...
            'CompareIndex',[1 2 3;1 2 3],'CompareColor',{condition_colors(1,:),condition_colors(2,:),condition_colors(3,:)},...
            'Colors',{condition_colors(1,:),condition_colors(2,:),condition_colors(3,:)},...
            'Transparent',transparent,'LineStyles',{'-'},'axes',h_1463_PSTH_diff, 'PCritical',0.01);
        
        axis tight; legend off;
        %         xlim([-300 2300]);
        %         xticks([0 1000 2050]); xticklabels({'Stim','1000 ms','Sac'});
        
        % Gaussian vel
        plot(xlim,[0 0], 'k--')
        % ylim([max(-20,min(ylim)) max(ylim)]);
        
        posi_neg_ratio = max(5, max(ylim)/abs(min(ylim)));
        %         ylim([-max(ylim)/posi_neg_ratio max(ylim)]);
        
        plot(Gauss_vel(:,1) + time_markers{1}(1), 0 + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
        
        title(sprintf('\\DeltaPSTH, p = %.3g,%.3g,%.3g',...
            group_result(ori_cell_no).ChoicePreference_pvalue(1,:)'));   % change from 3 to 1, ZZ 20210506
        
        % --- 3. Choice Divergence with newly added p value ---
        h_1463_CD = h_1463_axes(3);
        axes(h_1463_CD)
        
        for j = 1:3
            CD_this{j} = ChoiceDiv_All{j}(ori_cell_no,:,:);
            CD_perm_this{j} = 1.96 * ChoiceDiv_All_perm{j}.std(ori_cell_no,:,:);
            ps_CD_this{j} = squeeze(ChoiceDiv_All_perm{j}.p(ori_cell_no,:,:));
        end
        
        ylim([-0.2 1.1]);
        
        SeriesComparison({CD_this{1}, CD_this{2} CD_this{3}},...
            {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
            'Colors',{condition_colors(1,:),condition_colors(2,:),condition_colors(3,:)},...
            'OverridePs',{ps_CD_this{1}, ps_CD_this{2} ps_CD_this{3}},'ErrorBar',4,...
            'CompareIndex',[1 2 3 ; 1 2 3],'CompareColor',{condition_colors(1,:),condition_colors(2,:),condition_colors(3,:)},...
            'Transparent',transparent,'LineStyles',{'-'},'axes',h_1463_CD,'YLim', [-1.1 1.1], 'PCritical',0.01);
        
        legend off; plot(xlim,[0 0],'k--');
        
        title('CD');
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
        %         xticks([0 1000 2050]); xticklabels({'Stim','1000 ms','Sac'});
        yticks([-1 -0.5 0 0.5 1]);  yticklabels(num2str([-1 -0.5 0 0.5 1],'%0.1f\n'));
        
        
        % -- PSTH_diff (Hard and Easy) -- HH20160905
        %         %{
        figure(1466); clf;
        set(gcf,'uni','norm','pos',[0.013       0.099       0.286       0.352]);
        
        for j = 1:3
            PSTH_diff_hardeasy{j} = PSTH_hard_easy_raw_cellBtB4Bk{j}(ori_cell_no,:,1:2:end,:)...
                - PSTH_hard_easy_raw_cellBtB4Bk{j}(ori_cell_no,:,2:2:end,:);
        end
        
        SeriesComparison({reshape(PSTH_diff_hardeasy{1},1,[],6) reshape(PSTH_diff_hardeasy{2},1,[],6) reshape(PSTH_diff_hardeasy{3},1,[],6)},...
            {rate_ts{1} rate_ts{2} rate_ts{3} time_markers},...
            'Colors',hsv2rgb([2/3 .4 1; 2/3 1 1; 0 .4 1; 0 1 1; 1/3 .4 1; 1/3 1 1]),...
            'LineStyles',{'-','-'},'figN',1466);
        title('PSTH\_diff, Hard and Easy');
        axis tight; legend off; plot(xlim,[0 0],'k--');
        
        % Gaussian vel
        plot(Gauss_vel(:,1) + time_markers{1}(1),min(ylim) + Gauss_vel(:,2)*range(ylim)/4,'--','linew',1.5,'color',[0.6 0.6 0.6]);
        %         xlim([-300 2300]);
        %}
        
        % -- Dora Tuning -- HH20170327 @ UNIGE
        if ~isempty(dora_tuning_mean_each_cell)
            
            ind_in_partial = sum(select_for_partial(1:ori_cell_no)); % From original_ind to index in select_partial
            unique_heading = group_result(representative_cell).unique_heading;
            
            set(figure(1465),'name','Dora tuning'); clf;
            set(gcf,'uni','norm','pos',[0.425       0.057       0.567       0.607]);
            to_plot_tt = [1 4 5 6];
            [h_sub,~] = tight_subplot(3,length(to_plot_tt),[0.05 0.02]);
            zero_index = find(unique_heading == 0);
            
            for tt = 1:length(to_plot_tt)
                for k = 1:3
                    axes(h_sub(k+(tt-1)*3));
                    
                    hh_index = {1:zero_index zero_index:length(unique_heading)};
                    cc_marker = {'<','>'};
                    
                    for cc = LEFT:RIGHT
                        for hh = LEFT:RIGHT
                            
                            h = scatter(unique_heading(hh_index{hh}),...
                                dora_tuning_mean_each_cell(ind_in_partial,to_plot_tt(tt),k,cc,hh_index{hh}),...
                                dora_tuning_n_each_cell(ind_in_partial,to_plot_tt(tt),k,cc,hh_index{hh})*20,...
                                condition_colors(k,:),[cc_marker{cc}]); hold on;
                            
                            if cc == hh
                                set(h,'markerfacecolor',condition_colors(k,:));
                            else
                                set(h,'markerfacecolor','none');
                            end
                        end
                        
                        plot(unique_heading,squeeze(dora_tuning_mean_each_cell(ind_in_partial,to_plot_tt(tt),k,cc,:)),'color',condition_colors(k,:),'LineWid',2)
                    end
                    
                    if k == 1
                        title(partial_corr_timewins{to_plot_tt(tt),2});
                    end
                    
                end
            end
        end

        
        Plot_memsac([], [], ori_cell_no, h_1463_axes);
    end

    function Plot_memsac(~, ~, ori_cell_no, h_1463_axes)    % Plot associated mem-sac traces @HH20150426
        
        
        h_1463_memPSTH = h_1463_axes(4);
        h_1463_memPolar = h_1463_axes(5);
        
        try
            
            if ~isempty(group_result(ori_cell_no).mat_raw_MemSac)
                
                % plot_time_range = [2 3 4 6]; % 2:end
                plot_time_range = [3]; % Memory
                
                resp_mean = group_result(ori_cell_no).mat_raw_MemSac.resp_mean(plot_time_range);
                resp_err = group_result(ori_cell_no).mat_raw_MemSac.resp_err(plot_time_range);
                p = group_result(ori_cell_no).mat_raw_MemSac.p(plot_time_range);
                vectSum = group_result(ori_cell_no).mat_raw_MemSac.vectSum(plot_time_range);
                temporal_Slice = group_result(ori_cell_no).mat_raw_MemSac.temporal_Slice(plot_time_range,:);
                align_offsets = group_result(ori_cell_no).mat_raw_MemSac.align_offsets;
                align_markers = group_result(ori_cell_no).mat_raw_MemSac.align_markers;  % Desired markers: target onset & target offset & sac onset
                unique_heading = group_result(ori_cell_no).mat_raw_MemSac.unique_heading;
                
                %             PREF_target_location = group_result(ori_cell_no).mat_raw_PSTH.PREF_target_location; % Note this would be [-90,270]! HH20160906
                %             PREF_target_location = mod(PREF_target_location,360); % Changed to [0,360] HH20160906
                
                % Memsac PSTH: align to VisOn and SacOn. HH20180609
                MemSac_interp_PSTH = group_result(ori_cell_no).mat_raw_MemSac.MemSac_interp_PSTH([1 3]);
                MemSac_interp_PSTH_sem = group_result(ori_cell_no).mat_raw_MemSac.MemSac_interp_PSTH_sem([1 3]);
                MemSac_interp_PSTH_ts = group_result(ori_cell_no).mat_raw_MemSac.t_centers([1 3]);
                
                mem_sac_align_offset = mean(align_offsets); % This is time of [VisON, VisON, Saccade]
                mem_sac_time_marker = {mem_sac_align_offset - mem_sac_align_offset(1), mem_sac_align_offset - mem_sac_align_offset(3)};
                mem_sac_border = [mem_sac_align_offset(2)-100 -300]; % for SeriesComparison
                mem_sac_lims = [-300 mem_sac_align_offset(2)+700];
                
                MemSac_interp_locations = group_result(ori_cell_no).MemSac_interp_locations;
                
                [~,pref] = min(abs(group_PREF_target_location(ori_cell_no) - MemSac_interp_locations));
                %                     [~,pref] = min(abs(group_PREF_target_location_notebook(ori_cell_no) - MemSac_interp_locations));
                pref = mod(pref-1,length(MemSac_interp_locations)-1)+1;
                null = mod(pref + (length(MemSac_interp_locations)-1)/2 -1, length(MemSac_interp_locations)-1)+1; % The opposite position
                
                % Polar plot
                [~, polarOrder] = sort(max(cell2mat(resp_mean'),[],2),1,'descend'); % The largest goes first in the polar plot
                
                axes(h_1463_memPolar);
                title_text = '';
                
                for sliceN = polarOrder(:)'
                    
                    if strcmp(temporal_Slice{sliceN,5},''); continue; end
                    
                    if p(sliceN) < 0.05
                        sigMarker = '-'; wid = 2;
                    else
                        sigMarker = '-'; wid = 1;
                    end
                    
                    polarwitherrorbar([unique_heading/180*pi; unique_heading(1)/180*pi], ...
                        [resp_mean{sliceN}'; resp_mean{sliceN}(1)], (p(sliceN) < 0.05 || 1 ) * [resp_err{sliceN}' ; resp_err{sliceN}(1)],...
                        [temporal_Slice{sliceN,5} sigMarker], wid, 1);  % Simple mode = 1
                    
                    %                     set(polar([unique_heading/180*pi; unique_heading(1)/180*pi], ...
                    %                         [resp_mean{sliceN}'; resp_mean{sliceN}(1)],[temporal_Slice{sliceN,5} sigMarker]),'linewidth',2);
                    %
                    hold on;
                    
                    h_p = polar([vectSum(sliceN) vectSum(sliceN)]/180*pi,...
                        [0  max(ylim)],[temporal_Slice{sliceN,5} '-']);
                    
                    
                    %                     h_p = polar([vectSum(sliceN) vectSum(sliceN)]/180*pi,...
                    %                                 [max(cell2mat(resp_mean))*0.9 max(cell2mat(resp_mean))*1.3],[temporal_Slice{sliceN,5} '-']);
                    
                    if p(sliceN) < 0.05
                        %         h_p = polar([vectSum(sliceN) vectSum(sliceN)]/180*pi,[0 vectAmp(sliceN)],[temporal_Slice{sliceN,5} '-']);
                        set(h_p,'LineWidth',3);
                    else
                        set(h_p,'LineWidth',1);
                    end
                    %
                    % axes(axis_typing); axis off;
                    %                     text(0,- sliceN * 0.15+1.1,sprintf('%s:  \\itp_ = %4.2g',temporal_Slice{sliceN,4},p(sliceN)),'color',temporal_Slice{sliceN,5},...
                    %                         'FontSize',14 * (p(sliceN) < 0.05) + 7 * (p(sliceN) >= 0.05 || isnan(p(sliceN))));
                    drawnow;
                    title_text = [title_text sprintf('%s:  \\itp_ = %4.2g\n',temporal_Slice{sliceN,4},p(sliceN))];
                end
                
                % Target position
                axes(h_1463_memPolar);
                radi = max(xlim);
                pref_theta = MemSac_interp_locations(pref)/180*pi;
                plot(radi*cos(pref_theta),radi*sin(pref_theta),'ok','markersize',5,'linew',2,'markerfacecol','k');
                
                null_theta = MemSac_interp_locations(null)/180*pi;
                plot(radi*cos(null_theta),radi*sin(null_theta),'ok','markersize',5,'linew',2);
                
                text(min(xlim),min(ylim),title_text);
                
                
                
                % (2) The two which are nearest to actual PREF location. @HH20150524
                if isempty(group_result(ori_cell_no).TwoPtMemSac_PSTH) % Only if there is not TwoPtMemSac
                    
                    
                    for jjj = 1:2
                        ys{jjj} = shiftdim(MemSac_interp_PSTH{jjj}(:,[pref null]),-1);
                        sems{jjj} = MemSac_interp_PSTH_sem{jjj}(:,[pref null]);
                    end
                    
                    axes(h_1463_memPSTH); hold on;
                    
                    SeriesComparison({ys{1},ys{2}},...
                        {MemSac_interp_PSTH_ts{1},MemSac_interp_PSTH_ts{2},mem_sac_time_marker},...
                        'OverrideError',{sems{1},sems{2}},...
                        'ErrorBar', 2, 'CompareIndex',[],...
                        'Colors',{'k'},'hold', 1,'Ylim',[0 max(max(max((ys{:}))))], ...
                        'Transparent',transparent,'Border', mem_sac_border, 'LineStyles',{'-','--'},'axes',h_1463_memPSTH);
                    
                    text(0,max(ylim)*.88,'ActualFixP\_interp','color','k');
                end
            end
            
            % (3) The two in 2pt-memsac. @HH20150524
            TwoPtMemSac_PSTH = group_result(ori_cell_no).TwoPtMemSac_PSTH;
            
            if ~isempty(TwoPtMemSac_PSTH)
                
                % Memsac PSTH 2pt: align to VisOn and SacOn. HH20180609
                MemSac_interp_PSTH_ts = group_result(ori_cell_no).mat_raw_2ptMemSac.t_centers([1 3]);
                
                align_offsets = group_result(ori_cell_no).mat_raw_MemSac.align_offsets;
                mem_sac_align_offset = mean(align_offsets); % This is time of [VisON, VisON, Saccade]
                mem_sac_time_marker = {mem_sac_align_offset - mem_sac_align_offset(1), mem_sac_align_offset - mem_sac_align_offset(3)};
                mem_sac_border = [mem_sac_align_offset(2)-100 -300]; % for SeriesComparison
                mem_sac_lims = [-300 mem_sac_align_offset(2)+700];
                
                plot_aligns = [1 3];
                
                for jjj = 1:2
                    MemSac_interp_PSTH{jjj} = reshape...
                        ([group_result(ori_cell_no).mat_raw_2ptMemSac.result_PSTH_anne_mean{plot_aligns(jjj),:}],...
                        length(MemSac_interp_PSTH_ts{jjj}),[]);
                    MemSac_interp_PSTH_sem{jjj} = reshape...
                        ([group_result(ori_cell_no).mat_raw_2ptMemSac.result_PSTH_anne_sem{plot_aligns(jjj),:}],...
                        length(MemSac_interp_PSTH_ts{jjj}),[]);
                end
                
                % If HD_pref is on the left, flip 2pt-MemSac to let PREF go first
                if 90 <= group_PREF_target_location(ori_cell_no) && group_PREF_target_location(ori_cell_no) <=270
                    for jjj = 1:2
                        MemSac_interp_PSTH{jjj} = fliplr(MemSac_interp_PSTH{jjj});
                        MemSac_interp_PSTH_sem{jjj} = fliplr(MemSac_interp_PSTH_sem{jjj});
                    end
                end
                
                for jjj = 1:2
                    ys{jjj} = shiftdim(MemSac_interp_PSTH{jjj},-1);
                    sems{jjj} = MemSac_interp_PSTH_sem{jjj};
                end
                
                axes(h_1463_memPSTH); hold on;
                
                SeriesComparison({ys{1},ys{2}},...
                    {MemSac_interp_PSTH_ts{1},MemSac_interp_PSTH_ts{2},mem_sac_time_marker},...
                    'OverrideError',{sems{1},sems{2}},...
                    'ErrorBar', 2, 'CompareIndex',[],...
                    'Colors',{'k'}, 'hold', 1, 'Ylim',[0 max(max(max((ys{:}))))],...
                    'Transparent',transparent,'Border', mem_sac_border, 'LineStyles',{'-','--'},'axes',h_1463_memPSTH);
                
                text(0,max(ylim)*.81,'2pt','color',hsv2rgb([0.8 1 0.6]));
                
            end
            
            % Annotating
            % xlabel('Time to saccade onset (ms)');
            ylabel('Firing rate (Hz)');
            xlim(mem_sac_lims)
            legend off;
            
            % title(sprintf('PREF: %g, %s',group_PREF_target_location(ori_cell_no),group_result(ori_cell_no).cellID{2}{1}(30:end)));
            title(sprintf('%s',group_result(ori_cell_no).cellID{2}{1}(30:end)));
            
            for sliceN = 1:size(temporal_Slice,1)
                if temporal_Slice{sliceN,3} == 4 % Precise windows
                    plot([temporal_Slice{sliceN,1} temporal_Slice{sliceN,2}],[max(ylim) max(ylim)],temporal_Slice{sliceN,5},'linew',5);
                else  % Windows that is not so precise
                    % Mean shift
                    meanShift = mean(align_offsets(:,align_markers == temporal_Slice{sliceN,3})-align_offsets(:,align_markers==4),1);
                    plot( [meanShift meanShift],ylim,'k--','LineWidth',1);
                    plot(meanShift + [temporal_Slice{sliceN,1} temporal_Slice{sliceN,2}],[max(ylim) max(ylim)],temporal_Slice{sliceN,5},'linew',5);
                end
            end
            
            % Plot target manually
            plot([0 mem_sac_time_marker{1}(2)],[max(ylim) max(ylim)],'r','linew',5);
            
            xticks([0 1000 sum(abs(mem_sac_border))+100]); xticklabels({'TargOn','1000 ms','Sac'});
            ylim([0 max(ylim)])
            
            SetFigure(13)
            
        catch err
            err
            % keyboard
        end
    end


end


function [corr, rt, estrength, block_id, condi] = ...
    load_data(data_set, subject_id, show_progress)
% [corr, rt, estrength, block_id, condi] = load_data(data_set, subject_id)
% returns the data of subject subject_id from the given data_set. data_set
% can take one of the following values:
% 'vest_pilot' - pilot with vestibular-only stimulus, Gaussian velocity
% profile
% corr is the vector of correct (1) or wrong (0) response, rt is the
% reaction time in seconds, estrength is the evidence strength. Block id is
% the block id string associated with the trial. condition is the condition
% id for each trial, given that the experiment has one

% process arguments
if nargin < 3 || ~show_progress
    show_progress = false;
end

% get directory of file
data_root = fileparts(mfilename('fullpath'));

% valid data_set?
data_sets = {'vest_pilot', 'vest_pilot2', 'vest_pilot3', 'vest_tri', ...
    'vis_pilot', 'vis_pilot2', 'vest_vis_gauss', 'vest_vis_gauss_ext', ...
    'vest_comb', 'monk_pilot', 'vest_vis_gauss_ZZ'};  % Added by ZZ @20231008
if sum(strcmp(data_set, data_sets)) == 0
    error('Unknown data set "%s"', data_set);
end

% valid subject id?
if strcmp(data_set, 'vest_vis_gauss_ext') && strcmp(subject_id, 'subj03_cut')
    subject_id = 'subj03';
    cut_rts = true;
else
    cut_rts = false;
end
subj_ids = subject_ids(data_set);
if nargin < 2
    if length(subj_ids) > 1
        error('Subject id required for data set "%s"', data_set);
    else
        subject_id = subj_ids{1};
    end
elseif sum(strcmp(subject_id, subj_ids)) ~= 1
   error('Unknown subject id "%s" for data set "%s"', subject_id, data_set);
end


% is there a summary file? - if yes, load instead
summary_file = [data_root filesep data_set filesep subject_id '_summary.mat'];

if isempty(dir(summary_file))
    raw = csvread('D:\Paper_rawdata\Code\Labtools\OptimalMultisensoryDecisionMakingwithRT-main\data\ALL_trial_info_m10.csv', 1,1);
    corr = raw(:,5) == 2;  % choice: right, 1; left, 0; 
    rt = raw(:,4);     % reaction time 
    estrength = raw(:,3);   % headings
    block_id = raw(:,1);  % block index
    condi = raw(:,2);   % stim type
    
    save(summary_file, 'corr', 'rt', 'estrength', 'block_id', 'condi'); 
end

if ~isempty(dir(summary_file))
    if show_progress, fprintf('Loading pre-processed summary file\n'); end
    d = load(summary_file);
    corr = d.corr;
    rt = d.rt;
    estrength = d.estrength;
    block_id = d.block_id;
    condi = d.condi;
    if cut_rts
        valid_trials = rt <= 2;
        corr = corr(valid_trials);
        rt = rt(valid_trials);
        estrength = estrength(valid_trials);
        block_id = block_id(valid_trials);
        condi = condi(valid_trials);
    end
    return
end

% get list of files (i.e. blocks) and potentially exclude some of them
file_ids = subject_ids(data_set, subject_id);
exclude_first_n = 0;
file_ids = sort(file_ids);
% exclude first n files
file_ids = file_ids((exclude_first_n + 1):end);


if strcmp(data_set, 'vest_pilot')
    % Pilot data, vestibular stimulus only, Gaussian velocity profile
    data_dir = [data_root filesep data_set filesep 'Lab_' subject_id];
    % re-convert to original filenames
    block_ids = file_ids;
    for f_idx = 1:length(file_ids)
        file_ids{f_idx} = [data_dir filesep subject_id '_70_' ...
            file_ids{f_idx}(8:11) '_' file_ids{f_idx}(1:6) '_vestib.mat'];
    end
    % load data
    [corr, rt, estrength, block_id, condi] = ...
        read_human_data_files(file_ids, block_ids, 0, show_progress);
    
elseif strcmp(data_set, 'vest_pilot2')
    % Pilot data, vestibular stimulus only, Gaussian velocity profile
    data_dir = [data_root filesep data_set];
    % re-convert to original filenames
    block_ids = file_ids;
    for f_idx = 1:length(file_ids)
        file_ids{f_idx} = [data_dir filesep 'sub01_vestib_' ...
            file_ids{f_idx} '.mat'];
    end
    % load data
    [corr, rt, estrength, block_id, condi] = ...
        read_human_data_files(file_ids, block_ids, 0, show_progress);

elseif strcmp(data_set, 'vest_pilot3')
    % Pilot data, vestibular stimulus only, Gaussian velocity profile
    % for additional subject that will be added to 'vest_vis_gauss' later
    data_dir = [data_root filesep data_set filesep subject_id];
    % re-convert to original filenames
    block_ids = file_ids;
    for f_idx = 1:length(file_ids)
        file_ids{f_idx} = [data_dir filesep subject_id '_' ...
            file_ids{f_idx} '.mat'];
    end
    % load data
    [corr, rt, estrength, block_id, condi] = ...
        read_human_data_files(file_ids, block_ids, 0, show_progress);

elseif strcmp(data_set, 'vest_tri')
    % Vestibular only, triangular velocity profile, different accelerations
    % interleaved
    data_dir = [data_root filesep data_set filesep subject_id];
    % re-convert to original filenames
    block_ids = file_ids;
    if strcmp(subject_id, 'H'), subj_num = '831195';
    elseif strcmp(subject_id, 'J'), subj_num = '548214';
    else subj_num = '231139';
    end
    for f_idx = 1:length(file_ids)
        file_ids{f_idx} = [data_dir filesep subj_num '_' file_ids{f_idx} ...
            '_interleaved.mat'];
    end
    % load data
    [corr, rt, estrength, block_id, condi] = ...
        read_human_data_files(file_ids, block_ids, 1, show_progress);
    
elseif strcmp(data_set, 'vis_pilot')
    % visual only, gaussian velocity profile, different coherences
    % interleaved
    data_dir = [data_root filesep data_set filesep subject_id];
    % re-convert to original filenames
    block_ids = file_ids;
    for f_idx = 1:length(file_ids)
        file_ids{f_idx} = [data_dir filesep subject_id '_visual_' ...
            file_ids{f_idx} '.mat'];
    end
    % load data
    [corr, rt, estrength, block_id, condi] = ...
        read_human_data_files(file_ids, block_ids, 2, show_progress);

elseif strcmp(data_set, 'vis_pilot2')
    % visual only, gaussian velocity profile, different coherences
    % interleaved, subjects from 2nd/3rd batch
    data_dir = [data_root filesep data_set filesep subject_id];
    % re-convert to original filenames
    block_ids = file_ids;
    for f_idx = 1:length(file_ids)
        file_ids{f_idx} = [data_dir filesep subject_id '_visual_' ...
            file_ids{f_idx} '.mat'];
    end
    % load data
    [corr, rt, estrength, block_id, condi] = ...
        read_human_data_files(file_ids, block_ids, 2, show_progress);

elseif strcmp(data_set, 'vest_vis_gauss')
    % visual/vestibular, gaussian velocity profile, 3 different coherences
    data_dir = [data_root filesep data_set filesep subject_id];
    % re-convert to original filenames
    block_ids = file_ids;
    for f_idx = 1:length(file_ids)
        file_ids{f_idx} = [data_dir filesep subject_id '_interleaved_' ...
            file_ids{f_idx} '.mat'];
    end
    % load data
    [corr, rt, estrength, block_id, condi] = ...
        read_human_data_files(file_ids, block_ids, 3, show_progress);

elseif strcmp(data_set, 'vest_vis_gauss_ext')
    % visual/vestibular, gaussian velocity profile, 6 different coherences
    data_dir = [data_root filesep data_set filesep subject_id];
    % re-convert to original filenames
    block_ids = file_ids;
    for f_idx = 1:length(file_ids)
        file_ids{f_idx} = [data_dir filesep subject_id '_' ...
            file_ids{f_idx} '.mat'];
    end
    % load data
    [corr, rt, estrength, block_id, condi] = ...
        read_human_data_files(file_ids, block_ids, 4, show_progress);

elseif strcmp(data_set, 'vest_comb')
    % vestibular/combined conditions, different coherences
    data_dir = [data_root filesep data_set filesep subject_id];
    % re-convert to original filenames
    block_ids = file_ids;
    for f_idx = 1:length(file_ids)
        file_ids{f_idx} = [data_dir filesep subject_id '_vestibularcombined_' ...
            file_ids{f_idx} '.mat'];
    end
    % load data
    [corr, rt, estrength, block_id, condi] = ...
        read_human_data_files(file_ids, block_ids, 3, show_progress);
    
elseif strcmp(data_set, 'monk_pilot')
    % monkey data, vis/vestibular, gaussian 2s profile, single coherence
    data_dir = [data_root filesep data_set filesep subject_id];
    % re-convert to original filenames
    if strcmp(subject_id, 'mon01'), fname_start = 'm5c0r';
    else fname_start = 'm6c0r'; end
    block_ids = file_ids;
    for f_idx = 1:length(file_ids)
        file_ids{f_idx} = [data_dir filesep ...
            fname_start file_ids{f_idx} '.mat'];
    end
    % load data
    [corr, rt, estrength, block_id, condi] = ...
        read_monkey_data_files(file_ids, block_ids, 0, show_progress);
end

% make sure that estrength does not contain duplicate (very similar) es
ess = unique(estrength);
while length(ess) > 1
    sim_idx = find(abs(ess(2:end)-ess(1))<1e-10,1);
    if ~isempty(sim_idx)
        % ess(1) and ess(sim_idx+1) similar enough to unify
        estrength(estrength==ess(sim_idx+1))=ess(1);
        ess = [ess(1:sim_idx) ess((sim_idx+2):end)];
    else
        ess = ess(2:end);
    end
end

% create summary file for later use
if show_progress, fprintf('Writing summary file\n'); end
save(summary_file, 'corr', 'rt', 'estrength', 'block_id', 'condi');

% cut RT's, if requested
if cut_rts
    valid_trials = rt <= 2;
    corr = corr(valid_trials);
    rt = rt(valid_trials);
    estrength = estrength(valid_trials);
    block_id = block_id(valid_trials);
    condi = condi(valid_trials);
end


% ------------------------------------------------------------------------
function [corr, rt, estrength, block_id, condi] = ...
    read_human_data_files(file_names, block_ids, file_type, show_progress)
% file_type:
% 0 - just different headings
% 1 - different headings and different acceleration magnitutde
% 2 - different headings and different visual coherence
% 3 - different headings, vis only/vest only/vis-vest, 3 different vis coh
% 4 - different headings, vis only/vest only/vis-vest, 6 different vis coh

% loads the data in files of given list
files = length(file_names);
if show_progress, fprintf('%d Files to load\n', files); end

% allocate space
file_corr = cell(1, files);
file_nullt = cell(1, files);
file_rt = cell(1, files);
file_es = cell(1, files);
file_condi = cell(1, files);

% process files, one by one
for f_idx = 1:files
    file_name = file_names{f_idx};
    if show_progress
        fprintf('%02d: Processing %s...\n', f_idx, file_name);
    end
    % load data and get file statistics
    d = load(file_name);
    blocks = length(d.SavedInfo.Resp);
    trials_per_block = length(d.SavedInfo.Resp(1,1).corr);
    trials_per_file = blocks * trials_per_block;
    if show_progress, fprintf('%d blocks, %d trials per block\n', ...
            blocks, trials_per_block); end
    % store data in vectors, one per file
    file_corr{f_idx} = zeros(1, trials_per_file);
    file_nullt{f_idx} = zeros(1, trials_per_file);
    file_rt{f_idx} = zeros(1, trials_per_file);
    file_es{f_idx} = zeros(1, trials_per_file);
    file_condi{f_idx} = zeros(1, trials_per_file);
    file_trials = 0;
    for block_idx = 1:blocks
        % skip blocks that do not have full amount of trials
        if length(d.SavedInfo.Resp(block_idx).corr) ~= trials_per_block
            if show_progress, fprintf('Ignoring block %02d  - only %d trials\n', ...
                block_idx, length(d.SavedInfo.Resp(block_idx).corr)); end
            continue;
        end
        
        % process block
        from_idx = file_trials + 1;
        to_idx = from_idx + trials_per_block - 1;
        file_trials = file_trials + trials_per_block;
        
        file_corr{f_idx}(from_idx:to_idx) = d.SavedInfo.Resp(block_idx).corr;
        file_nullt{f_idx}(from_idx:to_idx) = d.SavedInfo.Resp(block_idx).null;
        if isfield(d.SavedInfo.Resp(block_idx), 'responeTime')
            % misspelling in first batch of data
            file_rt{f_idx}(from_idx:to_idx) = ...
                d.SavedInfo.Resp(block_idx).responeTime ...
                - d.SavedInfo.Resp(block_idx).delayTime;
        else
            file_rt{f_idx}(from_idx:to_idx) = ...
                d.SavedInfo.Resp(block_idx).responseTime ...
                - d.SavedInfo.Resp(block_idx).delayTime;            
        end
        file_es{f_idx}(from_idx:to_idx) = d.SavedInfo.Resp(block_idx).dir;
        
        % condition number depends on filetype
        if file_type == 1
            % acceleration value determines condition
            for trial = 1:trials_per_block
                accel = d.SavedInfo.Rep(block_idx).Trial(trial).Param(1).value;
                if abs(accel - 0.02) < 1e-10, file_condi{f_idx}(from_idx + trial - 1) = 1;
                elseif abs(accel - 0.04) < 1e-10, file_condi{f_idx}(from_idx + trial - 1) = 2;
                elseif abs(accel - 0.08) < 1e-10, file_condi{f_idx}(from_idx + trial - 1) = 3;
                elseif abs(accel - 0.1) < 1e-10, file_condi{f_idx}(from_idx + trial - 1) = 4;
                else
                    error('Unexpected acc value "%f", block %d, trial %d', ...
                        accel, block_idx, trial);
                end
            end
        elseif file_type == 2
            % coherence value determines condition
            for trial = 1:trials_per_block
                coh = d.SavedInfo.Rep(block_idx).Trial(trial).Param(30).value;
                if abs(coh - 25) < 1e-10, file_condi{f_idx}(from_idx + trial - 1) = 1;
                elseif abs(coh - 37) < 1e-10, file_condi{f_idx}(from_idx + trial - 1) = 2;
                elseif abs(coh - 70) < 1e-10, file_condi{f_idx}(from_idx + trial - 1) = 3;
                else
                    error('Unexpected coh value "%f", block %d, trial %d', ...
                        coh, block_idx, trial);
                end
            end
        elseif file_type == 3
            % on one hand coherence values, on the other hand coherences
            for trial = 1:trials_per_block
                % stimulus type encoded in condition as
                % 0x00 - 1 = visual stimulus present
                % x000 - 1 = vestibular stimulus presentplot_vis_pilot.m
                % i.e. 4 - vis, 8 - vest, 12 vis/vest
                s_type = d.SavedInfo.Rep(block_idx).Trial(trial).Param(52).value;
                if abs(s_type) < 1e-10, condi = 4;
                elseif abs(s_type - 1) < 1e-10, condi = 8;
                elseif abs(s_type - 2) < 1e-10, condi = 12;
                else
                    error('Unexpected stimulus type value "%f", block %d, trial %d', ...
                        s_type, block_idx, trial);
                end
                if bitand(condi, 4) > 0
                    % only care about coherence if visual stimulus present
                    % visual stimulus encoded in bits 00xx
                    % 1 - 25%, 2 - 37%, 3 - 70%
                    coh = d.SavedInfo.Rep(block_idx).Trial(trial).Param(30).value;
                    if abs(coh - 25) < 1e-10, condi = condi + 1;
                    elseif abs(coh - 37) < 1e-10, condi = condi + 2;
                    elseif abs(coh - 70) < 1e-10, condi = condi + 3;
                    else
                        error('Unexpected coh value "%f", block %d, trial %d', ...
                            coh, block_idx, trial);
                    end
                end
                file_condi{f_idx}(from_idx + trial - 1) = condi;
            end
        elseif file_type == 4
            % on one hand coherence values, on the other hand coherences
            for trial = 1:trials_per_block
                % stimulus type encoded in condition as
                % 0x000 - 1 = visual stimulus present
                % x0000 - 1 = vestibular stimulus presentplot_vis_pilot.m
                % i.e. 8 - vis, 16 - vest, 24 vis/vest
                s_type = d.SavedInfo.Rep(block_idx).Trial(trial).Param(52).value;
                if abs(s_type) < 1e-10, condi = bin2dec('01000');
                elseif abs(s_type - 1) < 1e-10, condi = bin2dec('10000');
                elseif abs(s_type - 2) < 1e-10, condi = bin2dec('11000');
                else
                    error('Unexpected stimulus type value "%f", block %d, trial %d', ...
                        s_type, block_idx, trial);
                end
                if bitand(condi, bin2dec('01000')) > 0
                    % only care about coherence if visual stimulus present
                    % visual stimulus encoded in bits 00xxx
                    % 1 - 0.1% (& 0%), 2 - 12%, 3 - 25%, 4 - 37%, 5 - 52%, 6 - 70%
                    coh = d.SavedInfo.Rep(block_idx).Trial(trial).Param(30).value;
                    if abs(coh) < 1e-10 || abs(coh - 0.1) < 1e-10, condi = condi + 1;
                    elseif abs(coh - 12) < 1e-10, condi = condi + 2;
                    elseif abs(coh - 25) < 1e-10, condi = condi + 3;
                    elseif abs(coh - 37) < 1e-10, condi = condi + 4;
                    elseif abs(coh - 52) < 1e-10, condi = condi + 5;
                    elseif abs(coh - 70) < 1e-10, condi = condi + 6;
                    else
                        error('Unexpected coh value "%f", block %d, trial %d', ...
                            coh, block_idx, trial);
                    end
                end
                file_condi{f_idx}(from_idx + trial - 1) = condi;
            end
        end
    end
    
    % remove ignored blocks
    file_corr{f_idx} = file_corr{f_idx}(1:file_trials);
    file_nullt{f_idx} = file_nullt{f_idx}(1:file_trials);
    file_rt{f_idx} = file_rt{f_idx}(1:file_trials);
    file_es{f_idx} = file_es{f_idx}(1:file_trials);
    file_condi{f_idx} = file_condi{f_idx}(1:file_trials);
    if show_progress, fprintf('%d trials loaded\n', file_trials); end
end

% pack the data of all files into single vector
corr = cell2mat(file_corr);
nullt = cell2mat(file_nullt);
rt = cell2mat(file_rt);
estrength = cell2mat(file_es);
condi = cell2mat(file_condi);
% create array that contains file_idx for each trials
block_num = cell2mat(cellfun(@(c) ones(1, length(file_corr{c})) * c, ...
    mat2cell(1:files, 1, ones(1, files)), 'UniformOutput', false));

% remove null trials
corr = corr(~logical(nullt));
rt = rt(~logical(nullt));
estrength = estrength(~logical(nullt));
block_num = block_num(~logical(nullt));
condi = condi(~logical(nullt));
% scale to s
rt = rt * 1e-3;
% create list of block id's
block_id = cell(size(block_num));
for n = 1:length(block_id);
    block_id{n} = block_ids{block_num(n)};
end


% ------------------------------------------------------------------------
function [corr, rt, estrength, block_id, condi] = ...
    read_monkey_data_files(file_names, block_ids, file_type, show_progress)
% file_type:
% 0 - just different headings, fixed coherence

% loads the data in files of given list
files = length(file_names);
if show_progress, fprintf('%d Files to load\n', files); end

% allocate space
file_corr = cell(1, files);
file_rt = cell(1, files);
file_es = cell(1, files);
file_condi = cell(1, files);

% process files one by one
for f_idx = 1:files
    file_name = file_names{f_idx};
    if show_progress
        fprintf('%02d: Processing %s...\n', f_idx, file_name);
    end
    % load data and get file statistics
    d = load(file_name);
    correct_trials = (d.rightward_trials & d.heading > 0) | ...
        (~d.rightward_trials & d.heading < 0);
    file_corr{f_idx} = double(correct_trials);
    file_rt{f_idx} = d.stim_durations;
    file_es{f_idx} = d.heading;
    file_condi{f_idx} = d.stim_type;
    if show_progress, fprintf('%d trials loaded\n', length(d.heading)); end
end

% pack the data of all files into single vector
corr = cell2mat(file_corr);
rt = cell2mat(file_rt);
estrength = cell2mat(file_es);
condi = cell2mat(file_condi);
% create array that contains file_idx for each trials
block_num = cell2mat(cellfun(@(c) ones(1, length(file_corr{c})) * c, ...
    mat2cell(1:files, 1, ones(1, files)), 'UniformOutput', false));

% scale to s
rt = rt * 1e-3;
% create list of block id's - turn sequential number into date
first_id = min(cellfun(@(c) str2double(c), block_ids));
for id_idx = 1:length(block_ids)
    % datestr seems to have a bug when handling 'mmddyy' format string -
    % work around that
    block_str = datestr(str2double(block_ids{id_idx}) - first_id, 'mm dd yy');
    block_ids{id_idx} = [block_str(1:2) block_str(4:5) block_str(7:8)];
end
block_id = cell(size(block_num));
for n = 1:length(block_id);
    block_id{n} = block_ids{block_num(n)};
end

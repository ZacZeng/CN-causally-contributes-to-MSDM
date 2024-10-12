function d = subject_ids(data_set, subject_id)
% d = subject_ids(data_set, subject_id)
% returns either a list of subject id's (if subject_id) is not given or
% subcategories for the given subject id and data set.

just_ids = (nargin == 1);
% data sets and subjects within data set
data_sets = {'vest_pilot', 'vest_pilot2', 'vest_pilot3', 'vest_tri', ...
    'vis_pilot', 'vis_pilot2', ...
    'vest_vis_gauss', 'vest_vis_gauss_ext', 'vest_comb', 'monk_pilot',...
    'vest_vis_gauss_ZZ'};
subj_ids = {{'J', 'M'}, {'J'}, {'950107', '967716'}, {'J', 'H', 'K'}, ...
    {'231139', '548203', '548214', '985010'}, ...
    {'950107', '967716', '541251', '950089'}, ...
    {'231139', '548203', '548214', '950107', '967716', ...
     '541251', '950089', '950091', '231140', '950096', '950129'}, ...
    {'subj01', 'subj02', 'subj03'}, ...
    {'950107'}, ...
    {'mon01', 'mon02'},...
    {'Monkey_M', 'RNN'}};   % Added by ZZ @20231008

% identify data set
data_set_id = find(strcmp(data_set, data_sets), 1);
if isempty(data_set_id)
    error('Unknown data set "%s"', data_set);
end

% identify subject ids
ids = subj_ids{data_set_id};
if just_ids
    d = ids;
    return
elseif sum(strcmp(subject_id, ids)) ~= 1
    error('Unknown subject id "%s" for data set "%s"', subject_id, data_set);
end

% get directory of file
data_root = fileparts(mfilename('fullpath'));

if strcmp(data_set, 'vest_pilot')
    % get file listing for given subject
    data_dir = [data_root filesep data_set filesep 'Lab_' subject_id];
    data_files = dir([data_dir filesep '*vestib.mat']);
    d = cell(1, length(data_files));
    % subid is 'mmddyy_runx'
    for d_idx = 1:length(d)
        d{d_idx} = ...
            [data_files(d_idx).name(11:16) '_' data_files(d_idx).name(6:9)];
    end
    
elseif strcmp(data_set, 'vest_pilot2')
    % get file listing for only subject
    data_dir = [data_root filesep data_set];
    data_files = dir([data_dir filesep 'subj01*.mat']);
    d = cell(1, length(data_files));
    % subid is 'mmddyy_runx'
    for d_idx = 1:length(d)
        d{d_idx} = data_files(d_idx).name(14:24);
    end

elseif strcmp(data_set, 'vest_pilot3')
    % get file listing for only subject
    data_dir = [data_root filesep data_set filesep subject_id];
    data_files = dir([data_dir filesep subject_id '_*.mat']);
    d = cell(1, length(data_files));
    % subid is 'mmddyy*' where * is everything before .mat
    for d_idx = 1:length(d)
        d{d_idx} = data_files(d_idx).name(...
            8:(findstr(data_files(d_idx).name, '.mat') - 1));
    end

elseif strcmp(data_set, 'vest_tri')
    % get file listing for given subject
    data_dir = [data_root filesep data_set filesep subject_id];
    data_files = dir([data_dir filesep '*interleaved.mat']);
    d = cell(1, length(data_files));
    % subid is 'mmddyy_runx'
    for d_idx = 1:length(d)
        d{d_idx} = data_files(d_idx).name(8:18);
    end
    
elseif strcmp(data_set, 'vis_pilot')
    % get file listing for given subject
    data_dir = [data_root filesep data_set filesep subject_id];
    data_files = dir([data_dir filesep '*_visual_*.mat']);
    d = cell(1, length(data_files));
    % subid is 'mmddyy_repx'
    for d_idx = 1:length(d)
        d{d_idx} = data_files(d_idx).name(15:25);
    end

elseif strcmp(data_set, 'vis_pilot2')
    % get file listing for given subject
    data_dir = [data_root filesep data_set filesep subject_id];
    data_files = dir([data_dir filesep '*_visual_*.mat']);
    d = cell(1, length(data_files));
    % subid is 'mmddyy_repx'
    for d_idx = 1:length(d)
        d{d_idx} = data_files(d_idx).name(15:25);
    end

elseif strcmp(data_set, 'vest_vis_gauss')
    % get file listing for given subject
    data_dir = [data_root filesep data_set filesep subject_id];
    data_files = dir([data_dir filesep '*_interleaved_*.mat']);
    d = cell(1, length(data_files));
    % subid is 'mmddyy_repx'
    for d_idx = 1:length(d)
        d{d_idx} = data_files(d_idx).name(20:30);
    end
    
elseif strcmp(data_set, 'vest_vis_gauss_ext')
    % get file listing for given subject
    data_dir = [data_root filesep data_set filesep subject_id];
    data_files = dir([data_dir filesep subject_id '_*_rep*.mat']);
    d = cell(1, length(data_files));
    % subid is 'mmddyy_hh.mm_repx'
        low_idx = length(subject_id) + 2;
    for d_idx = 1:length(d)
        d{d_idx} = data_files(d_idx).name(low_idx:(end-4));
    end
    
elseif strcmp(data_set, 'vest_comb')
    % get file listing for given subject
    data_dir = [data_root filesep data_set filesep subject_id];
    data_files = dir([data_dir filesep '*_vestibularcombined_*']);
    d = cell(1, length(data_files));
    % subid is 'mmddyy_repx'
    for d_idx = 1:length(d)
        d{d_idx} = data_files(d_idx).name(27:37);
    end
    
elseif strcmp(data_set, 'monk_pilot')
    % get file listing for given subject
    data_dir = [data_root filesep data_set filesep subject_id];
    data_files = dir([data_dir filesep 'm*c0r*.mat']);
    d = cell(1, length(data_files));
    % subid is 'xxxx' or 'xxxxx' - running number
    for d_idx = 1:length(d)
        d{d_idx} = data_files(d_idx).name(...
            6:(findstr(data_files(d_idx).name, '.') - 1));
    end
    
elseif strcmp(data_set, 'vest_vis_gauss_ZZ')   % Added by ZZ @ 20231008
    data_dir = [data_root filesep data_set filesep subject_id];
    data_files = dir([data_dir filesep 'm*c0r*.mat']);
    d = cell(1, length(data_files));
    % subid is 'xxxx' or 'xxxxx' - running number
    for d_idx = 1:length(d)
        d{d_idx} = data_files(d_idx).name(...
            6:(findstr(data_files(d_idx).name, '.') - 1));
    end
    
end
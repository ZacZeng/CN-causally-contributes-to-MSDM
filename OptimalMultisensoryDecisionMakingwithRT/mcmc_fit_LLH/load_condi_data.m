function d = load_condi_data(subj_id, show_progress)
%% loads all conditions of the data for a given subjects
%
% the data is returned blocked by condition
% vis - Nx visual
% vest - 1x vestibular
% comb - Nx combined
% where N is the number of visual coherences.
%
% if subj_id is a 6-digit number, then data is loaded from vest_vis_guass.
% Otherwise it is loaded from vest_vis_gauss_ext.

%% default argument
if nargin < 1
    subj_id = 'Monkey_M';
    show_progress = false;
end
if nargin < 2
    show_progress = false;
end
% % check if subject is from vest_vis_gauss
% if isempty(regexp(subj_id, '^\d{6}$', 'once'))
%     data_set = 'vest_vis_gauss_ext';
%     cohs = [0.001 0.12 0.25 0.37 0.52 0.70];
%     vis_base = 8;
%     vest_base = 16;
%     comb_base = 24;
% else
%     data_set = 'vest_vis_gauss';
%     cohs = [0.25 0.37 0.70];
%     vis_base = 4;
%     vest_base = 8;
%     comb_base = 12;
% end

% Changed by ZZ for my dataset
data_set = 'vest_vis_gauss_ZZ';
cohs = 0.3;
vis_base = 1;
vest_base = 1;
comb_base = 2;

coh_num = length(cohs);


%% load data
[corr, rt, es, block_id, condi] = load_data(data_set, subj_id, show_progress);
% for some subjects, remove some days
if strcmp(subj_id, 'subj03')
    % remove first 4 days of kalpana dataset
    [corr, rt, es, block_id, condi] = ...
        remove_days(corr, rt, es, block_id, condi, 4);
end


%% block data by conditions
% visual data
vis(1, coh_num) = struct('rt', [], 'choice', [], 'es', []);
for condi_id = 1:coh_num
    condi_trials = condi == vis_base + condi_id;
    condi_corr = corr(condi_trials);
    condi_rt = rt(condi_trials);
    condi_es = es(condi_trials);
    choice = logical(condi_corr); 
%     choice = (condi_es > 0 & logical(condi_corr)) | ...
%         (condi_es < 0 & ~logical(condi_corr));
    vis(condi_id) = struct('rt', condi_rt, 'choice', choice, 'es', condi_es);
end

% vestibular data
condi_trials = condi == vest_base;
condi_corr = corr(condi_trials);
condi_rt = rt(condi_trials);
condi_es = es(condi_trials);
choice = logical(condi_corr);
% choice = (condi_es > 0 & logical(condi_corr)) | (condi_es < 0 & ~logical(condi_corr));
vest = struct('rt', condi_rt, 'choice', choice, 'es', condi_es);

% combined data
comb(1, coh_num) = struct('rt', [], 'choice', [], 'es', []);
for condi_id = 1:coh_num
    condi_trials = condi == (condi_id + comb_base);
    condi_corr = corr(condi_trials);
    condi_rt = rt(condi_trials);
    condi_es = es(condi_trials);
    choice = logical(condi_corr); 
%     choice = (condi_es > 0 & logical(condi_corr)) | ...
%         (condi_es < 0 & ~logical(condi_corr));
    comb(condi_id) = struct('rt', condi_rt, 'choice', choice, 'es', condi_es);
end


%% combine into single structure
d = struct('vis', vis, 'vest', vest, 'comb', comb, 'cohs', cohs);



function [corr, rt, es, block_id, condi] = ...
    remove_days(corr, rt, es, block_id, condi, days_thresh)
%% returns the data with the first days_thresh days removed

% remove the first few days of trials
day0 = sort(unique(block_id));
day0 = datenum([2013 str2double(day0{1}(1:2)) str2double(day0{1}(3:4))]);
% index array picking only trials for which time is at least day_thresh
% away from first day
block_trials = arrayfun(@(i) (datenum([2013 ...
    str2double(block_id{i}(1:2)) str2double(block_id{i}(3:4))]) - day0) ...
    >= days_thresh, 1:length(block_id));

corr = corr(block_trials);
rt = rt(block_trials);
es = es(block_trials);
block_id = block_id(block_trials);
condi = condi(block_trials);

function [f_samp, m_samp] = family_wise(...
    log_ev, fam_model_idx, verbose, sample_num)
%% performs family-wise model comparison under RFX assumption
%
% If verbose (default false) is given and true, then a count of the samples
% is shown during the sampling process.
%
% sample_num defaults to 10000, if not given.
%
% The function returns the family probability samples, and the model
% probability samples, with samples being the rows of the returned
% matrices.


%% settings
if nargin < 4
    sample_num = 10000;
    if nargin < 3
        verbose = false;
    end
end
burnin = 10000;


%% process data
subj_num = size(log_ev, 1);
fam_num = length(fam_model_idx);
fam_model_num = arrayfun(@(i) length(fam_model_idx{i}), 1:fam_num);
model_num = sum(fam_model_num);


%% model prior parameters
model_pri = zeros(1, model_num);
for fam_idx = 1:fam_num
    for model_idx = 1:fam_model_num(fam_idx)
        % uniform distribution over each family
        model_pri(fam_model_idx{fam_idx}(model_idx)) = ...
            1/fam_model_num(fam_idx);
    end
end


%% initial samples
m = gamrnd(model_pri, 1);
m = m / sum(m);


%% re-sample
if verbose
    fprintf('Taking %d samples (+ %d burnin samples)\n', sample_num, burnin);
end
m_samp = zeros(burnin+sample_num, model_num);
for s_idx = 1:(burnin+sample_num)
    % count samples
    if verbose
        if mod(s_idx, 1000) == 0, fprintf('%7d ', s_idx); end
        if mod(s_idx, 10000) == 0, fprintf('\n'); end
    end
    
    % compute model posterior per subject, combining prior and evidence
    model_post = bsxfun(@plus, log_ev, log(m));
    model_post = exp(bsxfun(@minus, model_post, max(model_post, [], 2)));
    model_post = bsxfun(@rdivide, model_post, sum(model_post, 2));
    
    % resample subject <-> model assignment
    model_count = zeros(1, model_num);
    for subj_idx = 1:subj_num
        z = randpick(1, model_post(subj_idx,:));
        model_count(z) = model_count(z) + 1;
    end
    
    % resample model probability
    m = gamrnd(model_count + model_pri, 1);
    m = m / sum(m);
    m_samp(s_idx, :) = m;
end
if verbose, fprintf('\n'); end
m_samp = m_samp((burnin+1):end,:);


%% gather samples by family
f_samp = zeros(sample_num, fam_num);
for fam_idx = 1:fam_num
    for model_idx = 1:fam_model_num(fam_idx)
        f_samp(:,fam_idx) = f_samp(:,fam_idx) + ...
            m_samp(:,fam_model_idx{fam_idx}(model_idx));
    end
end

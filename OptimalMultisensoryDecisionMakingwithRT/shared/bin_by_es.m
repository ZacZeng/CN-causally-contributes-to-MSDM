function r = bin_by_es(choice, rt, es, varargin)
%% r = bin_by_rt(choice, rs, es, ...)
%
% bins reaction time and choices by evidence strength and returns
% various statistics per bin.
%
% Optional parameters are:
%
% bin_by_es(..., 'bias', bias) adds the bias to the evidence strength when
% evaluating correct / incorrect choices. The evidence strength itself is
% not modified. By default the bias is 0.
%
% bin_by_es(..., 'weight', weight) adds an additional weight to each of the
% data points. weight needs to be of same size as choice, rt, es. By
% default, a weight of 1 is used for each observation.
%
% bin_by_es(..., 'pconf', pconf) evaluates the (pconf * 100)% confidence
% intervals for probability choice/correct.
%
% bin_by_es(..., 'escorr', false) does not correct the choice values for very
% small evidence strengths.
%
% The following stats are returned:
% es: array of unique evidence strengths
% pc: probability correct for each es
% pr: probability choice for each es
% rt: mean reaction time (standard error rt_se)
% rt_corr: mean reaction time, conditioned on correct (with rt_corr_se)
% rt_incorr: mean reaction time, conditioned on incorrect (with rt_incorr_se)
% w: weight per es, in terms of observed samples, sum(w) = 1
% w_corr: weight correct choices per es
% w_incorr: weight incorrect choices per es, sum(w_corr)+sum(w_incorr) = 1
% pc_ivar, rt_ivar, rt_corr_ivar, rt_incorr_ivar:
%     inverse variance estimates of estimators
%
% If pconf is given, the additional bounds are returned
% pc_l, pc_u: lower and upper confidence interval, probability correct
% pr_l, pl_u: lower and upper confidence interval, probability choice
%
% The inverse variances of the estimators are based on the following: For
% pc, the n observations are assumed to follow a binomial distribution,
% such that the variance of pc is pc (1 - pc) / n. In case of rt, the
% variances are just the squared standard errors.


%% process arguments
p = inputParser;
p.FunctionName = 'bin_by_es';
p.addRequired('choice', @isvector);
p.addRequired('rt', @isvector);
p.addRequired('es', @isvector);
p.addParamValue('bias', 0, @isscalar);
p.addParamValue('weight', [], @isvector);
p.addParamValue('pconf', [], @isscalar);
p.addParamValue('escorr', true, @isscalar);
p.parse(choice, rt, es, varargin{:});

trials = length(choice);
if length(rt) ~= trials || length(es) ~= trials
    error('choice, rt, and es need to be of equal length');
end
choice = choice(:);  rt = rt(:);  es = es(:);
bias = p.Results.bias;
if isempty(p.Results.weight)
    weight = ones(1, length(rt));
else
    weight = p.Results.weight(:)';
    if length(weight) ~= trials
        error('weight need to have same length as choice');
    end
end
pconf = p.Results.pconf;
escorr = p.Results.escorr;


%% process data
rt_raw = rt(:)';
es_raw = es(:)';
choice = choice(:)';
% get unique evidence strength
es = unique(es_raw);
es_num = length(es);
% allocate space
rt = zeros(1, es_num);             rt_se = zeros(1, es_num);
rt_corr = zeros(1, es_num);        rt_corr_se = zeros(1, es_num);
rt_incorr = zeros(1, es_num);      rt_incorr_se = zeros(1, es_num);
pr = zeros(1, es_num);             pc = zeros(1, es_num);
w_corr = zeros(1, es_num);         w_incorr = zeros(1, es_num);
w = zeros(1, es_num);
pc_ivar = zeros(1, es_num);         rt_ivar = zeros(1, es_num);
rt_corr_ivar = zeros(1, es_num);    rt_incorr_ivar = zeros(1, es_num);
if ~isempty(pconf)
    pr_l = zeros(1, es_num);   pr_u = zeros(1, es_num);
    pc_l = zeros(1, es_num);   pc_u = zeros(1, es_num);
end
% threshold on evidence strength discr = 0.1% of maximum absolute value
es_thresh = 0.001 * max(abs(es));
% correct choices - depends on bias
biased_es = es_raw + bias;
correct = biased_es >= es_thresh & logical(choice) | ...
    biased_es <= -es_thresh & ~logical(choice);
% evidence below threshold - correct 50% of time
es_thresh_idx = biased_es < es_thresh & biased_es > -es_thresh;
if escorr
    correct(es_thresh_idx) = ...
        0.5 * ones(1, sum(es_thresh_idx)) <= rand(1, sum(es_thresh_idx));
else
    correct(es_thresh_idx) = logical(choice(es_thresh_idx));
end

% bin data
for es_idx = 1:es_num
    d_idx = (es_raw == es(es_idx));
    d_corr_idx = d_idx & correct;
    d_incorr_idx = d_idx & ~correct;
    es_w = weight(d_idx);               w(es_idx) = sum(es_w);
    es_w_corr = weight(d_corr_idx);     w_corr(es_idx) = sum(es_w_corr);
    es_w_incorr = weight(d_incorr_idx); w_incorr(es_idx) = sum(es_w_incorr);
    % probability right & correct
    pr(es_idx) = sum(choice(d_idx) .* es_w) ./ w(es_idx); % NaN if empty
    pc(es_idx) = w_corr(es_idx) / w(es_idx);
    pc_ivar(es_idx) = w(es_idx) / (pc(es_idx) * (1 - pc(es_idx)));
    if ~isempty(pconf)
        if any(choice(d_idx))
            [~, aa] = ...
                binofit(sum(choice(d_idx) .* es_w), round(w(es_idx)), pconf);
            pr_l(es_idx)=aa(1); pr_u(es_idx)=aa(2); 
            [~, bb] = ...
                binofit(w_corr(es_idx), round(w(es_idx)), pconf);
            pc_l(es_idx)=bb(1); pc_u(es_idx)=bb(2); 
%             [pr_l(es_idx) pr_u(es_idx)] = ...
%                 binoconf(sum(choice(d_idx) .* es_w), round(w(es_idx)), pconf);
%             [pc_l(es_idx) pc_u(es_idx)] = ...
%                 binoconf(w_corr(es_idx), round(w(es_idx)), pconf);
        else
            pr_l(es_idx) = NaN;    pr_u(es_idx) = NaN;
            pc_l(es_idx) = NaN;    pc_u(es_idx) = NaN;
        end
    end
    % reaction times
    rt(es_idx) = sum(es_w .* rt_raw(d_idx)) / w(es_idx);
    rt_se(es_idx) = sqrt(sum(es_w .* ...
        (rt_raw(d_idx) - rt(es_idx)) .^ 2)) / w(es_idx);
    rt_ivar(es_idx) = 1 / rt_se(es_idx) ^ 2;
    rt_corr(es_idx) = sum(es_w_corr .* rt_raw(d_corr_idx)) / w_corr(es_idx);
    rt_corr_se(es_idx) = sqrt(sum(es_w_corr .* ...
        (rt_raw(d_corr_idx) - rt_corr(es_idx)) .^ 2)) / w_corr(es_idx);
    rt_corr_ivar(es_idx) = 1 / rt_corr_se(es_idx) ^ 2;
    rt_incorr(es_idx) = sum(es_w_incorr .* rt_raw(d_incorr_idx)) / w_incorr(es_idx);
    rt_incorr_se(es_idx) = sqrt(sum(es_w_incorr  .* ...
        (rt_raw(d_incorr_idx) - rt_incorr(es_idx)) .^ 2)) / w_incorr(es_idx);
    rt_incorr_ivar(es_idx) = 1 / rt_incorr_se(es_idx) ^ 2;
end

% fix NaN from 0 right/correct choices
pr(isnan(pr)) = 0;
pc(isnan(pc)) = 0;
pc_ivar(isnan(pc_ivar)) = 0;    pc_ivar(isinf(pc_ivar)) = 0;
% fix NaN for reaction times, from no available samples
rt_corr(isnan(rt_corr)) = 0;        rt_corr_se(isnan(rt_corr_se)) = 0;
rt_incorr(isnan(rt_corr)) = 0;      rt_incorr_se(isnan(rt_incorr_se)) = 0;
rt_corr_ivar(isnan(rt_corr_ivar)) = 0;
rt_corr_ivar(isinf(rt_corr_ivar)) = 0;
rt_incorr_ivar(isnan(rt_incorr_ivar)) = 0;
rt_incorr_ivar(isinf(rt_incorr_ivar)) = 0;

% normalise weights
n = w;
w_corr = w_corr ./ sum(w);
w_incorr = w_incorr ./ sum(w);
w = w ./ sum(w);

% return arrays
r = struct('es', es, 'pc', pc, 'pr', pr, 'pc_ivar', pc_ivar, ...
    'rt', rt, 'rt_se', rt_se, 'rt_ivar', rt_ivar, ....
    'rt_corr', rt_corr, 'rt_corr_se', rt_corr_se, 'rt_corr_ivar', rt_corr_ivar, ... 
    'rt_incorr', rt_incorr, 'rt_incorr_se', rt_incorr_se, 'rt_incorr_ivar', rt_incorr_ivar, ...
    'w', w, 'w_corr', w_corr, 'w_incorr', w_incorr, 'n', n);
if ~isempty(pconf)
    r.pr_l = pr_l;  r.pr_u = pr_u;
    r.pc_l = pc_l;  r.pc_u = pc_u;
end
function fit_model(subj_id, model_name, sample_num, ini_type)
%% fits the model to data for the given subject
%
% If ini_type is given, it can be one of
%
% 'resample' - only valid if previous data exists. In this case, it will
%     load the best-found llh of the previous dataset and use it as the
%     starting point, but will discard the previous samples otherwise.
%
% 'nograd' - only valid if no previous data exists. In this case, the
%     function will NOT perform an initial gradient ascent procedure before
%     sampling, but instead start sampling from the standard initialisation
%     point.


%% detect if running on cluster, process arguments accordingly
on_cluster = ~usejava('desktop');
if on_cluster
    fprintf('Desktop not detected - assuming running on cluster\n');
    sample_num = str2double(sample_num);
    % needs to be adjusted if you want to use this on a cluster
    base_folder = '/home/jdrugowitsch/vis_vest/mcmc_fit_LLH';
else
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======= Change the file directory===========
    base_folder = 'D:\Paper_rawdata\Code\Labtools\OptimalMultisensoryDecisionMakingwithRT-main\data\vest_vis_gauss_ZZ';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
if nargin >= 4 && ~isempty(ini_type)
    if strcmp(ini_type, 'resample')
        ini_type = 1;
    elseif strcmp(ini_type, 'nograd')
        ini_type = 2;
    else
        error('Unknown ini_type %s', ini_type);
    end
else ini_type = 0;
end
fprintf('Fitting model %s, subject %s, %d samples\n', ...
    model_name, subj_id, sample_num);


%% settings
data_folder = [base_folder filesep 'fit_data'];
model_fn = str2func(model_name);
data_filename = [data_folder filesep model_name '_' subj_id '.mat'];
opt=optimset('Display','iter','MaxFunEval',10000,'FunValCheck','on', ...
    'MaxIter',5000,'Algorithm','interior-point','TolX',1e-20);
delta_t = 5e-3;
t_max = 1.5;
% t_max = 2;


%% load subject data
if on_cluster
    d = load([base_folder filesep 'subj_data' filesep subj_id '.mat']);
else
    d = load_condi_data(subj_id, true);
end
dstats = get_data_stats(d);
n = dstats.n;
% add additional stats to d to be used by the model
d.hs = dstats.vis(1).hs;
d.delta_t = delta_t;
d.t_max = t_max;
d.subj_id = subj_id;
% ts = (0:ceil(t_max / delta_t))' * delta_t;
% [d.v_t, ~, d.a_t,~] = fit_real_vel_profile(delta_t);   d.a_t = abs(1000*d.a_t); 
[raw_vel, d.v_t, raw_acc, d.a_t] = fit_real_vel_profile(delta_t);   
% d.v_t = raw_vel; 
% d.a_t = raw_acc*1000; 
d.a_t = abs(d.a_t); 

%% get model parameter bounds
pstruct = model_fn(d);
p_ini = pstruct.p_ini;
p_min = pstruct.p_min;
p_max = pstruct.p_max;
p_w = pstruct.p_w;
p_names = pstruct.p_names;
p_num = length(p_ini);


%% log-likelihood function to minimize
logpdf = @(p) model_LLH(p, model_fn, d);


%% load previous fit data, if exists
s = zeros(sample_num,p_num);
llhs = zeros(sample_num,1);
s_base_idx = 0;
if exist(data_filename, 'file')
    % found previous data file
    if ini_type == 2
        error('ini_type nograd cannot be used if previous data exists');
    end
    fprintf('Loading previous fit data from %s\n', data_filename);
    m = load(data_filename);
    if ini_type == 1
        % discard previous samples, but use previous starting point
        p_ini = m.best_p;
        fprintf('Initialised with previous best parameters, but discarded samples\n');
        fprintf('Initial LLH=%f\n', m.best_llh);
    else
        % add to previous samples
        s = [m.s; zeros(sample_num,p_num)];
        llhs = [m.llhs; zeros(sample_num,1)];
        s_base_idx = size(m.s,1);
        p_ini = m.s(end,:);
        fprintf('Loaded %d previous samples\n', s_base_idx);
    end
else
    % no previous data found
    if ini_type == 1
        error('ini_type resample cannot be used without previous data');
    end
    fprintf('No previous fit data found\n');
    if ini_type == 2
        fprintf('Starting sampling from initial point\n');
    else
        fprintf('Initial gradient descent to find sampling starting point\n');
        [p_ini, llh] = fmincon(@(p) bounded_f(...
            @(x) -logpdf(x), p, p_min, p_max, 'exp'), ...
            p_ini, [], [], [], [], p_min, p_max, [], opt);
        fprintf('Best initial point found at LLH=%f\n', -llh);
    end
end


%% initialise sampler
sls = slsample_init(logpdf, p_ini, p_min, p_max, p_w);


%% take samples according to model function
fprintf('\nTaking %d samples...\n', sample_num);
s_tic = tic;
for s_idx = 1:sample_num
    % output
    if mod(s_idx, 10) == 0, fprintf('%d ', s_idx); end
    if mod(s_idx, 100) == 0
        s_time = toc(s_tic);
        fprintf('%7.3f samples/s\n', 100/s_time);
        s_tic = tic;
    end
    
    % sampling
    [s(s_idx+s_base_idx,:), sls] = slsample(sls, 1);
    llhs(s_idx+s_base_idx) = sls.logpdf_z;
end

[~,s_best] = max(llhs);
fprintf('\n\nFound best fit after sampling with LLH=%f at\n', llhs(s_best));
for p_idx = 1:length(p_names)
    fprintf('%s = %f\n', p_names{p_idx}, s(s_best,p_idx));
end
fprintf('\n');

    
%figure;
%for p_idx = 1:param_num
%    subplot(param_num+1,1,p_idx);
%    plot(samples(:,p_idx));
%    ylabel(param_names{p_idx});
%    set(gca,'Layer','top','TickDir','out','XTick',[],'XTickLabel',{});
%end
%subplot(param_num+1,1,param_num+1)
%plot(samples(:,end));
%ylabel('llh');
%set(gca,'Layer','top','TickDir','out');


%% from the best sample, perform another gradient ascent
[~,s_best] = max(llhs);
fprintf('\n\nPerforming gradient ascent from best sample, LLH=%f\n', llhs(s_best));
[p, llh] = fmincon(@(p) bounded_f(...
        @(x) -logpdf(x), p, p_min, p_max, 'exp'), ...
        s(s_best,:), [], [], [], [], p_min, p_max, [], opt);
best_llh = -llh;

fprintf('Found best fit with LLH=%f at\n', best_llh);
for p_idx = 1:length(p_names)
    fprintf('%s = %f\n', p_names{p_idx}, p(p_idx));
end
best_p = p;


%% re-computing model fits at this point
fprintf('\n\nRecomputing model fits\n');
[mstats, biases] = model_fn(d, best_p);
dstats = get_data_stats(d, biases);


%% store data
cohs = d.cohs;
fprintf('Writing data to %s\n', data_filename);
save(data_filename, 's', 'llhs', 'cohs', ...
    'best_llh', 'best_p', 'p_names', 'model_name', 'dstats', 'mstats', ...
    'subj_id');


%% exit matlab if running on cluster
if on_cluster
    exit;
end


function l = model_LLH(p, model_fn, d)
%% returns the LLH for the given parameters, model and subject data
[ms, biases] = model_fn(d, p);
ds = get_data_stats(d, biases);
l = get_LLH(ds, ms);


function r2 = model_R2(p, model_fn, d)
%% returns the R2 for the given parameters, model and subject data
[ms, biases] = model_fn(d, p);
ds = get_data_stats(d, biases);
[rt_var, pr_var] = get_R2_vars(ds);
[rt_err, pr_err] = get_R2_errs(ds, ms);
r2 = 1 - 0.5 * (rt_err/rt_var + pr_err/pr_var);
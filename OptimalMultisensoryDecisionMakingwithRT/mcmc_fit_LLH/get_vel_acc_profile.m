function [v_t, a_t] = get_vel_acc_profile(delta_t)
%% returns the velocity and acceleration profile

%% profile parameters
% t_max = 2; % path duration
% s_max = 0.3; % distance travelled
% num_sig = 4; % parameter settings of real velocity profile
% sigv = t_max / (8^(1/4) * num_sig); % conversion to actualy standard deviation
t_max = 1.5; % path duration
s_max = 0.12; % distance travelled
num_sig = 4.5; % parameter settings of real velocity profile
sigv = t_max / (8^(1/4) * num_sig); % conversion to actualy standard deviation


%% generate profiles
ts = (0:ceil(t_max / delta_t))' * delta_t;
t0 = t_max / 2;
s_t = sqrt(2 * pi * sigv^2) * ...
    (normcdf((ts - t0) / sigv) - normcdf(-t0 / sigv));
% adjust maximum velocity (i.e. profile gain) s.t. s_max correct
v_max = s_max / s_t(end);
%s_t = s_t * v_max;
v_t = v_max * exp(- 0.5 / sigv^2 * (ts - t0).^2);
a_t = ((t0 - ts) / sigv^2) .* v_t;

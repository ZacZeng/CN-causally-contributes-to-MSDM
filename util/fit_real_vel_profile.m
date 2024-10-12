% This script fits the normalized real speed and acceleration profiles for further analysis
function [raw_vel, nor_vel, raw_acc, nor_acc] = fit_real_vel_profile(delta_t)

% delta_t: time steps of time period, such as 0.005s
% raw_vel, raw_acc: real velocity and acceleration 
% nor: function of normalized velocity and abs(acceleration)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======= Change the file directory===========

% Load real velocity profile
load('D:\Paper_rawdata\Raw_data\CN\Gaussian_vel_real_sigma045amp12.mat')
% Time
t = Gaussian_vel_real_sigma045amp12(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    delta_t = t; 
elseif ~isempty(delta_t)
    wt = 1500;  % whole stimulus period = 1.5s
    delta_t = linspace(0, wt, ceil(1.5 / delta_t));
end

% Normalizing velocity profile
vel = Gaussian_vel_real_sigma045amp12(:,2);
nor_vel = vel / max(vel); 

%% Fitting velocity profile with Gaussian profile
gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
startPoints = [1 922 20 0.1]; 

f_vel = fit(t, vel, gaussEqn, 'Start', startPoints);  % Real velocity
f_nor_vel = fit(t, nor_vel, gaussEqn, 'Start', startPoints);   % Normilized velocity

%% Derivative of the above function as Acceleration 
d_f = ((2* f_vel.a *f_vel.b - 2* f_vel.a * t)/(f_vel.c ^2 )) .* exp(-((t-f_vel.b)/f_vel.c) .^2);
d_f_nor = ((2* f_nor_vel.a *f_nor_vel.b - 2* f_nor_vel.a * t)/(f_nor_vel.c ^2 )) .* exp(-((t-f_nor_vel.b)/f_nor_vel.c) .^2); 

% % Nomalization 
abs_nor_acc = d_f_nor / max(abs(d_f_nor));

% abs_nor_accEqn = 'abs(((a - b*x)*exp(-((x-c)/d)^2))/e)'; 
nor_accEqn = '((a - b*x)*exp(-((x-c)/d)^2))/e'; 
acc_startPoints = [0.023 0.00003 913 276 0.003];

f_acc = fit(t, d_f,nor_accEqn, 'Start', acc_startPoints);
f_nor_acc = fit(t, abs_nor_acc, nor_accEqn, 'Start', acc_startPoints);


%% Plotting
if nargin < 1
    
    set(figure(818), 'name', 'Normalized Velocity and Acceleration Profile', 'pos', [10 10 1500 500]); clf;
    
    % Vel
    subplot(1,2,1)
    plot(f_nor_vel, t, nor_vel)
    xlabel('Time (ms)'); ylabel('Normalized Velocity');
    
    % Acc
    subplot(1,2,2)
    plot(f_nor_acc, t, d_f_nor)
    xlabel('Time (ms)'); ylabel('Normalized Acceleration');
    
end
%% 
raw_vel = f_vel(delta_t); 
raw_acc = f_acc(delta_t);

nor_vel = f_nor_vel(delta_t);
nor_acc = f_nor_acc(delta_t);
return; 



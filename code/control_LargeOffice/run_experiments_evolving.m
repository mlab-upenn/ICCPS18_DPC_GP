% This script runs all the tracking experiment simulations to get data.
% It uses evolving GP method.

% clear all;

simindex = 1;

% Battery parameters
Batt_params = struct('power_max', 100000, ...% Battery max charge/discharge power [W]
    'power_min', -100000, ...   % Battery min charge/discharge power [W]
    'soc_max', 100000, ...     % Battery max SoC [Wh]
    'soc_min', 0, ...       % Battery min SoC [Wh]
    'power_weight', 0.85, ...   % probability of satisfying the power constraint
    'soc_weight', 0.85, ...     % probability of satisfying the SoC constraint
    'timestep', 0.25,...    % time step of the model
    'soc_discount', 1);  % discount rate of the weight soc_weight

running_time_stats = [];
reftol = 0.1; % 0.05*abs(DRreduction);

ramplimit = 2;

evolving_params = struct('period', 96*7,...   % number of steps before updating the model
    'method', 'IG', ...     % method: IG = Information Gain (only method)
    'maxpoints', 2016);      % max number of points allowed

StartMonth = 8;
StartDay = 4;  % Monday


%% Select the GP model
%control_gp_model_file = 'random_uniform_ramped2_ahead00_21days_truongkernel';
control_gp_model_file = 'random_uniform_ramped2_ahead00_21days_truongisokernel';
%control_gp_model_file = 'random_uniform_ramped2_ahead00_7days_truongkernel';

%%

DRstart = 14;
DRend = 16;
DRreduction = 90*1000;  % in Watts
%reftol = 0.05*abs(DRreduction);

% Set to true if multiple GPs are used for predictions; otherwise a single
% GP is used iteratively.
multipleGPs = false;

horizon = 4;
active_ahead = 4;

use_battery = true;
wdelta = 100;
wPb = 5;
%wvar = 0.2*wdelta;
wvar = 0.2;

use_evolving = true;

MATpostfix = sprintf('_battery_evolving_%s', datestr(now, 'yyyymmdd_HHMM'));  % _%02d  simindex

t = tic;
run_sim;
running_time_stats(simindex) = toc(t);
delete(controller)
clear controller
simindex = simindex + 1;


%% Summary
disp('Running time statistics:');
disp(running_time_stats);

% This script runs all the tracking experiment simulations to get data.

% clear all;

simindex = 1;

% Battery parameters
Batt_params = struct('power_max', 100000, ...% Battery max charge/discharge power [W]
    'power_min', -100000, ...   % Battery min charge/discharge power [W]
    'soc_max', 50000, ...     % Battery max SoC [Wh]
    'soc_min', 0, ...       % Battery min SoC [Wh]
    'power_weight', 0.9, ...   % weight for the variance in the power constraint
    'soc_weight', 0.95, ...     % weight or minimum weight for the variance in the SoC constraint
    'soc_weight_max', 1.6449, ... % max weight for the variance in the SoC constraint
    'timestep', 0.25,...    % time step of the model
    'soc_discount', 1);  % discount rate of the weight soc_weight

running_time_stats = [];
reftol = 0.1; % 0.05*abs(DRreduction);

% Values for weights above:
%   Probability of satisfying       Weight
%   p                               erfinv(p)*sqrt(2)
%   0.3                             0.3853
%   0.4                             0.5244
%   0.5                             0.6745
%   0.6                             0.8416
%   0.7                             1.0364
%   0.75                            1.1503
%   0.8                             1.2816
%   0.9                             1.6449
%   0.95                            1.96 ~ 2

%% Select the GP model
control_gp_model_file = 'doe_noreset_IG_2ramped_3input_ahead00_3days_truongkernel';
%control_gp_model_file = 'random_uniform_2ramped_3input_ahead00_3days_truongkernel';

ramplimit = 2;

%{
%%

DRstart = 13;
DRend = 16;
DRreduction = 30*1000;  % in Watts
%reftol = 0.05*abs(DRreduction);

% Set to true if multiple GPs are used for predictions; otherwise a single
% GP is used iteratively.
multipleGPs = false;

horizon = 8;
active_ahead = 4;

use_battery = true;
wdelta = 100;
wPb = 0;
%wvar = 0.2*wdelta;
wvar = 10;

MATpostfix = sprintf('_%s', datestr(now, 'yyyymmdd_HHMM'));  % _%02d  simindex

t = tic;
run_sim;
running_time_stats(simindex) = toc(t);
delete(controller)
clear controller
simindex = simindex + 1;
%}

%%

DRstart = 13;
DRend = 16;
DRreduction = 20*1000;  % in Watts
%reftol = 0.05*abs(DRreduction);

% Set to true if multiple GPs are used for predictions; otherwise a single
% GP is used iteratively.
multipleGPs = false;

horizon = 4;
active_ahead = 4;

use_battery = false;
wdelta = 0;
wPb = 100;
%wvar = 0.2*wdelta;
wvar = 0;

MATpostfix = sprintf('_nobattery_%s', datestr(now, 'yyyymmdd_HHMM'));  % _%02d  simindex

t = tic;
run_sim;
running_time_stats(simindex) = toc(t);
delete(controller)
clear controller
simindex = simindex + 1;


%% Summary
disp('Running time statistics:');
disp(running_time_stats);

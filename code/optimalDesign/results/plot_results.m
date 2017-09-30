
no_inputs = 3;
no_days = 7;
const = 'constrained_'; % or ""

%% OED information gain reset
% load(['doe_sampling_reset_IG_' num2str(no_inputs) 'input_' num2str(no_days) 'day.mat'])
% LP(LP<-5) = NaN;
% 
% figure(1);  hold on; grid on; box on;
% h1igr = plot(LP_map, 'LineWidth', 2);
% ylabel('log probability')
% 
% figure(2);  hold on; grid on; box on;
% h2igr = plot(RMSE_map, 'LineWidth', 2);
% ylabel('RMSE')

%% OED information gain noreset
% load(['doe_sampling_noreset_' const 'IG_' num2str(no_inputs) 'input_' num2str(no_days) 'day.mat'])
load doe_sampling_noreset_constrained_IG_3input_15day_20170921_2330

LP(LP<-5) = NaN;

figure(1);  hold on; grid on; box on;
h1ig = plot(LP_map, 'LineWidth', 2);
ylabel('log probability')

figure(2);  hold on; grid on; box on;
h2ig = plot(RMSE_map, 'LineWidth', 2);
ylabel('RMSE')

%% RANDOM uniform

% load(['random_sampling_uniform_' const num2str(no_inputs) 'input_' num2str(no_days) 'day.mat'])
load random_sampling_uniform_constrained_3input_15day_20170921_2023

LP(LP<-5) = NaN;

figure(1)
h1un = plot(LP, '*', 'LineWidth', 2);
% legend([h1igr, h1ig, h1un], 'ig reset', 'ig noreset', 'uniform', 'Location', 'Best')
legend([h1ig, h1un], 'ig noreset', 'uniform', 'Location', 'Best')

figure(2)
h2un = plot(RMSE, '*', 'LineWidth', 2);
% legend([h2igr, h2ig, h2un], 'ig reset', 'ig noreset', 'uniform', 'Location', 'Best')
legend([h2ig, h2un], 'ig noreset', 'uniform', 'Location', 'Best')

%% RANDOM prbs

% load(['random_sampling_prbs_' num2str(no_inputs) 'input_' num2str(no_inputs) 'day.mat'])
% LP(LP<-5) = NaN;
% 
% figure(1)
% h1prbs = plot(LP, 'LineWidth', 2);
% legend([h1ig, h1un, h1prbs], 'ig', 'uniform', 'prbs', 'Location', 'Best')
% 
% figure(2)
% h2prbs = plot(RMSE, 'LineWidth', 2);
% legend([h2ig, h2un, h2prbs], 'ig', 'uniform', 'prbs', 'Location', 'Best')

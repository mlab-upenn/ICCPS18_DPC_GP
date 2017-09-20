
no_inputs = 1;
no_days = 10;

%% OED information gain reset
load(['doe_sampling_reset_IG_' num2str(no_inputs) 'input_' num2str(no_days) 'day.mat'])
LP(LP<-5) = NaN;

figure(1);  hold on; grid on; box on;
h1igr = plot(LP_map, 'LineWidth', 2);
ylabel('log probability')

figure(2);  hold on; grid on; box on;
h2igr = plot(RMSE_map, 'LineWidth', 2);
ylabel('RMSE')

%% OED information gain noreset
load(['doe_sampling_noreset_IG_' num2str(no_inputs) 'input_' num2str(no_days) 'day.mat'])
LP(LP<-5) = NaN;

figure(1);  hold on; grid on; box on;
h1ig = plot(LP_map, 'LineWidth', 2);
ylabel('log probability')

figure(2);  hold on; grid on; box on;
h2ig = plot(RMSE_map, 'LineWidth', 2);
ylabel('RMSE')

%% RANDOM uniform

load(['random_sampling_uniform_' num2str(no_inputs) 'input_' num2str(no_days) 'day.mat'])
LP(LP<-5) = NaN;

figure(1)
h1un = plot(LP, 'LineWidth', 2);
legend([h1igr, h1ig, h1un], 'ig reset', 'ig noreset', 'uniform', 'Location', 'Best')

figure(2)
h2un = plot(RMSE, 'LineWidth', 2);
legend([h2igr, h2ig, h2un], 'ig reset', 'ig noreset', 'uniform', 'Location', 'Best')

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

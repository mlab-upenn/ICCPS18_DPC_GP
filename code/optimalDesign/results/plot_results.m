

%% OED information gain
load doe_sampling_IG_1input_7day.mat
LP(LP<-5) = NaN;

figure(1);  hold on; grid on; box on;
h1ig = plot(LP, 'LineWidth', 2);
ylabel('log probability')

figure(2);  hold on; grid on; box on;
h2ig = plot(RMSE, 'LineWidth', 2);
ylabel('RMSE')

%% RANDOM uniform

load random_sampling_uniform_1input_7day.mat
LP(LP<-5) = NaN;

figure(1)
h1un = plot(LP, 'LineWidth', 2);

figure(2)
h2un = plot(RMSE, 'LineWidth', 2);

%% RANDOM prbs

load random_sampling_prbs_1input_7day.mat
LP(LP<-5) = NaN;

figure(1)
h1prbs = plot(LP, 'LineWidth', 2);
legend([h1ig, h1un, h1prbs], 'ig', 'uniform', 'prbs', 'Location', 'Best')

figure(2)
h2prbs = plot(RMSE, 'LineWidth', 2);
legend([h2ig, h2un, h2prbs], 'ig', 'uniform', 'prbs', 'Location', 'Best')
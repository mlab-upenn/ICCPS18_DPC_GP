
rng(0);

%% laod data

file = 'unconstrained-LargeHotel';
order_autoreg = 3;
n_samples = 1000;

ctrl_variables = {'ClgSP', 'KitchenClgSP', 'GuestClgSP', 'SupplyAirSP', 'ChwSP'};
[X, y] = load_data(file, order_autoreg, ctrl_variables);
X_train = X(1:n_samples,:);
y_train = y(1:n_samples);
X_test = X(n_samples+1:end,:);
y_test = y(n_samples+1:end);

% standardize the data set
[X_train_norm, X_train_min, X_train_max] = preNorm(X_train);
[y_train_norm, y_train_min, y_train_max] = preNorm(y_train);

X_test_norm = preNorm(X_test, X_train_min, X_train_max);
y_test_norm = preNorm(y_test, y_train_min, y_train_max);

%% define model

model.mean_function       = {@constant_mean};
model.covariance_function = {@ard_sqdexp_covariance};
model.likelihood          = @likGauss;

% used saved initial hyperparameters
load('init_hyp.mat');

% uncomment to calculate new initial hyperparams
n_samples_init = 1000;
init_hyp = initial_model(file, n_samples_init, order_autoreg, ctrl_variables);

true_hyp = init_hyp;

% priors on each log covariance parameter
priors.cov = cell(1,numel(init_hyp.cov));
for idc = 1:numel(init_hyp.cov)
    priors.cov{idc}  = get_prior(@gaussian_prior, init_hyp.cov(idc), 1);
%     priors.cov{idc}  = get_prior(@gaussian_prior, 0, 1);
end

% prior on log noise
priors.lik  = {get_prior(@gaussian_prior, init_hyp.lik, 1)};
% priors.lik  = {get_prior(@gaussian_prior, 0, 1)};

% prior on constant mean
priors.mean = {get_prior(@gaussian_prior, init_hyp.mean, 1)};
% priors.mean = {get_prior(@gaussian_prior, 0, 1)};

model.prior = get_prior(@independent_prior, priors);
model.inference_method = add_prior_to_inference_method(@exact_inference, model.prior);

num_points = size(X_train,1);
x_star = X_train_norm;
y_star = y_train_norm;

% setup problem struct
problem.num_evaluations  = 100;
problem.candidate_x_star = x_star;

% function is a simple lookup table
problem.f = @(x) (y_star(find(all(bsxfun(@eq, x, x_star), 2))));

% actively learn GP hyperparameters
results = learn_gp_hyperparameters(problem, model);

[f_star_mean_active, f_star_variance_active, ~, ~, log_probabilities] = ...
    gp(results.map_hyperparameters(end), model.inference_method, ...
       model.mean_function, model.covariance_function, model.likelihood, ...
       results.chosen_x, results.chosen_y, X_test_norm, y_test_norm);
f_star_mean_active = postNorm(f_star_mean_active, y_train_min, y_train_max);
f_star_variance_active = postNormVar(f_star_variance_active, y_train_min, y_train_max);

report_active = sprintf('ACTIVE:\n E[log p(y* | x*, D)] = %0.3f, RMSE = %0.1f', ...
                 mean(log_probabilities), sqrt(mean((f_star_mean_active-y_test).^2)));
fprintf('%s\n', report_active);

loss(y_test, f_star_mean_active, f_star_variance_active);

X_chosen_active = results.chosen_x;
y_chosen_active = results.chosen_y;
X_chosen_active = postNorm(X_chosen_active, X_train_min, X_train_max);
y_chosen_active = postNorm(y_chosen_active, y_train_min, y_train_max);

% compare to random sampling
ind = randperm(num_points, problem.num_evaluations);

X_chosen = x_star(ind, :);
y_chosen = y_star(ind);

map_hyperparameters_random = minimize_minFunc(model, X_chosen, y_chosen);

[f_star_mean, f_star_variance, ~, ~, log_probabilities] = ...
    gp(map_hyperparameters_random, model.inference_method, ...
       model.mean_function, model.covariance_function, model.likelihood, ...
       X_chosen, y_chosen, X_test_norm, y_test_norm);
f_star_mean = postNorm(f_star_mean, y_train_min, y_train_max);
f_star_variance = postNormVar(f_star_variance, y_train_min, y_train_max);


report = sprintf('RANDOM:\n E[log p(y* | x*, D)] = %0.3f, RMSE = %0.1f', ...
                 mean(log_probabilities), sqrt(mean((f_star_mean-y_test).^2)));
fprintf('%s\n', report);
loss(y_test, f_star_mean, f_star_variance);

X_chosen = postNorm(X_chosen, X_train_min, X_train_max);
y_chosen = postNorm(y_chosen, y_train_min, y_train_max);

% plotgp for active learning
t = [0:length(y_test)-1]';
f1=figure('Name', 'active learning');
f1 = plotgp(f1, t, y_test, f_star_mean_active, sqrt(f_star_variance_active));
axis1 = findobj(f1,'Type','axes');
axis1(2).XLim = [0 100];
axis1(1).XLim = [0 100];

% plotgp for random sampling
t = [0:length(y_test)-1]';
f2=figure('Name', 'random sampling');
f2 = plotgp(f2, t, y_test, f_star_mean, sqrt(f_star_variance));
axis1 = findobj(f2,'Type','axes');
axis1(2).XLim = [0 100];
axis1(1).XLim = [0 100];


% plot with input on x axis
idx = 14; % ChwSP

figure;
subplot(2, 1, 1);
hold on;
plot(X_train(:,idx), y_train, 'r.', 'MarkerSize', 5);
% plot(X_test(:,idx), f_star_mean_active, 'm.', 'MarkerSize', 5);
plot(X_chosen_active(:,idx), y_chosen_active, 'k+', 'MarkerSize', 15);
% axis([22, 32, 0, 5e5]);
title(report_active);

subplot(2, 1, 2);
hold on;
plot(X_train(:,idx), y_train, 'r.', 'MarkerSize', 5);
% plot(X_train(:,idx), f_star_mean, 'm.', 'MarkerSize', 5);
plot(X_chosen(:,idx), y_chosen, 'k+', 'MarkerSize', 15);
% axis([22, 32, 0, 5e5]);
title(report);
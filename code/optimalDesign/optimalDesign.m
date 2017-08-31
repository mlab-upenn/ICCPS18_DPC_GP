
rng(0);

%% laod data

file = 'unconstrained-LargeHotel';
ctrl_horizon = 1;
order_autoreg = 3;

[X, y] = load_data(file, order_autoreg, ctrl_horizon);
X_train = X;
y_train = y;

% standardize the data set
[X_train_norm, X_train_min, X_train_max] = preNorm(X_train);
[y_train_norm, y_train_min, y_train_max] = preNorm(y_train);


%% define model

model.mean_function       = {@constant_mean};
model.covariance_function = {@ard_sqdexp_covariance};
model.likelihood          = @likGauss;

% used saved initial hyperparameters
load('init_hyp.mat');
true_hyp = init_hyp;

% uncomment to calculate new initial hyperparams
% n_samples_init = 1000;
% init_hyp = initial_model(file, n_samples_init);

% priors on each log covariance parameter
priors.cov  = ...
    {get_prior(@gaussian_prior, 0, 1), ...
     get_prior(@gaussian_prior, 0, 1), ...
     get_prior(@gaussian_prior, 0, 1), ...
     get_prior(@gaussian_prior, 0, 1), ...
     get_prior(@gaussian_prior, 0, 1), ...
     get_prior(@gaussian_prior, 0, 1), ...
     get_prior(@gaussian_prior, 0, 1), ...
     get_prior(@gaussian_prior, 0, 1), ...
     get_prior(@gaussian_prior, 0, 1), ...
     get_prior(@gaussian_prior, 0, 1), ...
     get_prior(@gaussian_prior, 0, 1), ...
     get_prior(@gaussian_prior, 0, 1), ...
     get_prior(@gaussian_prior, 0, 1), ...
     get_prior(@gaussian_prior, 0, 1), ...
     get_prior(@gaussian_prior, 0, 1)};

% prior on log noise
priors.lik  = {get_prior(@gaussian_prior, 0, 1)};

% prior on constant mean
priors.mean = {get_prior(@gaussian_prior, 0, 1)};

model.prior = get_prior(@independent_prior, priors);
model.inference_method = ...
    add_prior_to_inference_method(@exact_inference, model.prior);

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

[~, ~, f_star_mean_active, f_star_variance_active, log_probabilities] = ...
    gp(results.map_hyperparameters(end), model.inference_method, ...
       model.mean_function, model.covariance_function, model.likelihood, ...
       results.chosen_x, results.chosen_y, x_star, y_star);
f_star_mean_active = postNorm(f_star_mean_active, y_train_min, y_train_max);
f_star_variance_active = postNormVar(f_star_variance_active, y_train_min, y_train_max);

report_active = sprintf('ACTIVE:\nE[log p(y* | x*, D)] = %0.3f, RMSE = %0.1f', ...
                 mean(log_probabilities), sqrt(mean((f_star_mean_active-y_train).^2)));
fprintf('%s\n', report_active);
loss(y_train, f_star_mean_active, f_star_variance_active);

X_chosen_active = results.chosen_x;
y_chosen_active = results.chosen_y;
X_chosen_active = postNorm(X_chosen_active, X_train_min, X_train_max);
y_chosen_active = postNorm(y_chosen_active, y_train_min, y_train_max);

% compare to random sampling
ind = randperm(num_points, problem.num_evaluations);

X_chosen = x_star(ind, :);
y_chosen = y_star(ind);

map_hyperparameters_random = minimize_minFunc(model, X_chosen, y_chosen);

[~, ~, f_star_mean, f_star_variance, log_probabilities] = ...
    gp(map_hyperparameters_random, model.inference_method, ...
       model.mean_function, model.covariance_function, model.likelihood, ...
       X_chosen, y_chosen, x_star, y_star);
f_star_mean = postNorm(f_star_mean, y_train_min, y_train_max);
f_star_variance = postNormVar(f_star_variance, y_train_min, y_train_max);


report = sprintf('RANDOM:\nE[log p(y* | x*, D)] = %0.3f, RMSE = %0.1f', ...
                 mean(log_probabilities), sqrt(mean((f_star_mean-y_train).^2)));
fprintf('%s\n', report);
loss(y_train, f_star_mean, f_star_variance);

X_chosen = postNorm(X_chosen, X_train_min, X_train_max);
y_chosen = postNorm(y_chosen, y_train_min, y_train_max);


idx = 14;

figure;
subplot(2, 1, 1);
hold on;
plot(X_train(:,idx), y_train, 'r.', 'MarkerSize', 5);
plot(X_train(:,idx), f_star_mean_active, 'm.', 'MarkerSize', 5);
plot(X_chosen_active(:,idx), y_chosen_active, 'k+', 'MarkerSize', 15);
% axis([22, 32, 0, 5e5]);
title(report_active);

subplot(2, 1, 2);
hold on;
plot(X_train(:,idx), y_train, 'r.', 'MarkerSize', 5);
plot(X_train(:,idx), f_star_mean, 'm.', 'MarkerSize', 5);
plot(X_chosen(:,idx), y_chosen, 'k+', 'MarkerSize', 15);
% axis([22, 32, 0, 5e5]);
title(report);


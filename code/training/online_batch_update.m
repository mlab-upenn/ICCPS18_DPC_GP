function sampled_data = online_batch_update(model, data)

% INPUTS:
% data should be struct with fields:
%             Ambient
%             Humidity
%             TotalLoad
%             TOD
%             DOW
%             ClgSP
%             KitchenClgSP
%             GuestClgSP
%             SupplyAirSP
%             ChwSP
%
% model should be struct with fields: 
%             init_hyp: initial hyperparams for OED
%             n_samples: batch size
%             model_inputs
%             model_target
%             model_stepsahead
% 
% OUTPUT:
% sampled_data is a struct with fields:
%             X: new features for training
%             y: new target values for training

%% load data

data_train = data;
normalized_fields = {'Ambient', 'Humidity', 'TotalLoad', 'ClgSP', 'KitchenClgSP', 'GuestClgSP', 'SupplyAirSP', 'ChwSP'};
[data_train_norm, normparams] = normalize_data(data_train, normalized_fields);

n_samples = model.n_samples;
model_inputs = model.model_inputs;
model_target = model.model_target;
stepsahead = model.model_stepsahead;
model_excepts = model.model_excepts;

[X_train_norm, y_train_norm] = construct_data(data_train_norm, model_inputs, model_target, stepsahead, model_excepts);

%% define model

model.mean_function       = {@constant_mean};
model.covariance_function = {@ard_sqdexp_covariance};
model.likelihood          = @likGauss;

% initial hyperparams
init_hyp = model.init_hyp;

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

x_star = X_train_norm;
y_star = y_train_norm;

% setup problem struct
problem.num_evaluations  = n_samples;
problem.candidate_x_star = x_star;

% function is a simple lookup table
problem.f = @(x) (y_star(find(all(bsxfun(@eq, x, x_star), 2))));

% actively learn GP hyperparameters
results = learn_gp_hyperparameters(problem, model);

% print error metrics on training data
print_results = 0;
if print_results
    [f_star_mean_active, f_star_variance_active, ~, ~, log_probabilities] = ...
        gp(results.map_hyperparameters(end), model.inference_method, ...
        model.mean_function, model.covariance_function, model.likelihood, ...
        results.chosen_x, results.chosen_y, X_train_norm, y_train_norm);
    f_star_mean_active = postNorm(f_star_mean_active, y_train_min, y_train_max);
    f_star_variance_active = postNormVar(f_star_variance_active, y_train_min, y_train_max);
    
    report_active = sprintf('ACTIVE:\n E[log p(y* | x*, D)] = %0.3f, RMSE = %0.1f', ...
        mean(log_probabilities), sqrt(mean((f_star_mean_active-y_test).^2)));
    fprintf('%s\n', report_active);
    loss(y_test, f_star_mean_active, f_star_variance_active);
end

X_chosen_active = results.chosen_x;
y_chosen_active = results.chosen_y;

% find the index of chosen samples in training data
ind = zeros(1,n_samples);
for id1 = 1:n_samples
    y_cur = y_chosen_active(id1);
    for id2 = 1:size(y_train_norm, 1)
        if y_cur == y_train_norm(id2)
            ind(id1) = id2;
            break
        end
    end
end

%% post process

% shift ind by however number of rows were removed from orginal data
ind = ind + size(data_train.Ambient,1) - size(X_train_norm,1);

sampled_data.indices = ind;
sampled_data.Ambient = data.Ambient(ind);
sampled_data.Humidity = data.Humidity(ind);
sampled_data.TotalLoad = data.TotalLoad(ind);
sampled_data.TOD = data.TOD(ind);
sampled_data.DOW = data.DOW(ind);
sampled_data.ClgSP = data.ClgSP(ind);
sampled_data.KitchenClgSP = data.KitchenClgSP(ind);
sampled_data.GuestClgSP = data.GuestClgSP(ind);
sampled_data.SupplyAirSP = data.SupplyAirSP(ind);
sampled_data.ChwSP = data.ChwSP(ind);


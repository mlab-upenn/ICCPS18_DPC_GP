%% load data

n_days_init = 14;
n_days_next = 14;
samples_init = 1:n_days_init*96;                        % June
shift = 30;                                             % no of daysto shift for next batch
samples_next = (shift)*96+1:(shift)*96+n_days_next*96;  % July
batch_size = 14;
n_samples = numel(samples_init)+numel(samples_next);

% slice data
building = 'LargeOffice';
datafname = 'random_LargeOffice_uniform_2ramped_3input_153day_20171005_1620.mat';
data_train = load(fullfile('data', datafname));
data_train = slice_data_batch(data_train, [samples_init, samples_next]);

% Normalize the data (not all fields)
normalized_fields = {'Ambient', 'Humidity', 'TotalLoad', 'ClgSP', 'LgtSP', 'SupplyAirSP', 'ChwSP'};
[data_train_norm, normparams] = normalize_data(data_train, normalized_fields);
y_train_min =  normparams.TotalLoad.min;
y_train_max = normparams.TotalLoad.max;

datafname = 'test-unconstrained-LargeOffice'; %'test-LargeOffice; 'test-ramped2-LargeOffice'
data_test = load(fullfile('..', 'data', datafname));

% Normalize the data (same as for training)
data_test_norm = normalize_data(data_test, normalized_fields, normparams);

model_inputs = {...
        'TOD', ...
        'DOW', ...
        {'Ambient', 1:-1:0}, ...
        {'Humidity', 1:-1:0}, ...
        {'TotalLoad', 1:-1:1}, ...
        {'ClgSP', 1:-1:0}, ... 
        {'LgtSP', 1:-1:0}, ...
        {'SupplyAirSP', 1:-1:0}, ...
        {'ChwSP', 1:-1:0}};
    
model_target = 'TotalLoad';
model_excepts = {'TOD', 'DOW'};
stepsahead = 0;

[X_train_norm, y_train_norm] = construct_data(data_train_norm, model_inputs, model_target, stepsahead, model_excepts);
[X_train, y_train] = construct_data(data_train, model_inputs, model_target, stepsahead, model_excepts);
[X_test_norm, y_test_norm] = construct_data(data_test_norm, model_inputs, model_target, stepsahead, model_excepts);
[~, y_test] = construct_data(data_test, model_inputs, model_target, stepsahead, model_excepts);

%% random sampling

ind = randperm(n_samples, batch_size);

X_chosen = X_train_norm(ind, :);
y_chosen = y_train_norm(ind);

model = train_gp(X_chosen, y_chosen);

[f_star_mean, f_star_variance, ~, ~, log_probabilities] = ...
    gp(model.hyp, model.inference_method, ...
       model.mean_function, model.covariance_function, model.likelihood, ...
       X_chosen, y_chosen, X_test_norm, y_test_norm);
f_star_mean = postNorm(f_star_mean, y_train_min, y_train_max);
f_star_variance = postNormVar(f_star_variance, y_train_min, y_train_max);

report = sprintf('RANDOM:\n E[log p(y* | x*, D)] = %0.3f, RMSE = %0.1f', ...
                 mean(log_probabilities), sqrt(mean((f_star_mean-y_test).^2)));
fprintf('%s\n', report);
loss(y_test, f_star_mean, f_star_variance);

X_chosen = X_train(ind, :);
y_chosen = y_train(ind);


clear;
clc;

rng(1);

%% load data

n_days_init = 21;
n_days_next = 21;
shift = 30; 
n_days = 16;
suffix = '';  % '_prior' or ''

loadStr = sprintf('evolving_%dinit_%dshift_%dupdate_%dday%s.mat',...
            n_days_init, shift, n_days_next, n_days, suffix);
load(fullfile('data',loadStr));


n_days_test = 14;
samples_init = 1:n_days_init*96;                            % June
                                                % no of days to shift for next batch
samples_next = (shift)*96+1:(shift)*96+n_days_next*96;
samples_test = samples_next(end) + (1:n_days_test*96);                    % test on next week of update batch
batch_size = n_days_init*96;
n_samples = numel(samples_init)+numel(samples_next);

savedata.n_days_init = n_days_init;
savedata.n_days_next = n_days_next;
savedata.shift = shift;

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

%% random sampling

new_batch_size = n_days*96;
ind = 1:new_batch_size;

X_chosen_random = X_train_norm(ind, :);
y_chosen_random = y_train_norm(ind);

model_random = train_gp(X_chosen_random, y_chosen_random);

%% active learning

model_active = train_gp(X_chosen_active, y_chosen_active);

%% print results for different data sets

data_test = load(fullfile('data', datafname));
data_test = slice_data_batch(data_test, samples_test);

% normalize the data (same as for training)
data_test_norm = normalize_data(data_test, normalized_fields, normparams);

[X_test_norm, y_test_norm] = construct_data(data_test_norm, model_inputs, model_target, stepsahead, model_excepts);
[~, y_test] = construct_data(data_test, model_inputs, model_target, stepsahead, model_excepts);

% random
[mu_random, s2_random, ~, ~, lp_random] = ...
    gp(model_random.hyp, model_random.inference_method, ...
    model_random.mean_function, model_random.covariance_function, model_random.likelihood, ...
    X_chosen_random, y_chosen_random, X_test_norm, y_test_norm);
mu_random = postNorm(mu_random, y_train_min, y_train_max);
s2_random = postNormVar(s2_random, y_train_min, y_train_max);

report = sprintf('RANDOM:\n E[log p(y* | x*, D)] = %0.3f, RMSE = %0.1f', ...
    mean(lp_random), sqrt(mean((mu_random-y_test).^2)));
fprintf('%s\n', report);
[aer, ser, lpdr, mrser, smser, msllr] = loss(y_test, mu_random, s2_random);

% plotgp
t = [0:length(y_test)-1]';
f=figure('Name', 'random sampling');
f = plotgp(f, t, y_test, mu_random, sqrt(s2_random));
axis1 = findobj(f,'Type','axes');
axis1(2).XLim = [0 size(X_test_norm,1)];
axis1(1).XLim = [0 size(X_test_norm,1)];

% active
[mu_active, s2_active, ~, ~, lp_active] = ...
    gp(model_active.hyp, model_active.inference_method, ...
    model_active.mean_function, model_active.covariance_function, model_active.likelihood, ...
    X_chosen_active, y_chosen_active, X_test_norm, y_test_norm);
mu_active = postNorm(mu_active, y_train_min, y_train_max);
s2_active = postNormVar(s2_active, y_train_min, y_train_max);

report_active = sprintf('ACTIVE:\n E[log p(y* | x*, D)] = %0.3f, RMSE = %0.1f', ...
    mean(lp_active), sqrt(mean((mu_active-y_test).^2)));
fprintf('%s\n', report_active);

[aea, sea, lpda, mrsea, smsea, mslla] = loss(y_test, mu_active, s2_active);

% plotgp
t = [0:length(y_test)-1]';
f=figure('Name', 'active learning');
f = plotgp(f, t, y_test, mu_active, sqrt(s2_active));
axis1 = findobj(f,'Type','axes');
axis1(2).XLim = [0 size(X_test_norm,1)];
axis1(1).XLim = [0 size(X_test_norm,1)];

%% save results data

% active
training_result = [];
training_result.stepsahead = stepsahead;
training_result.model_inputs = model_inputs;
training_result.model_target = model_target;
training_result.model_excepts = model_excepts;
training_result.hyp = model_active.hyp;
training_result.model = model_active;
Xtrain_norm = X_chosen_active;
Ytrain_norm = y_chosen_active;

validation_result = [];
validation_result.loss_ae = aea;
validation_result.loss_se = sea;
validation_result.loss_lpd = lpda;
validation_result.loss_mrse = mrsea;
validation_result.loss_smse = smsea;
validation_result.loss_msll = mslla;
validation_result.loss_rmse = sqrt(sea);

saveStr = sprintf('evolving_%dinit_%dshift_%dupdate_%dday%s.mat',...
            n_days_init, shift, n_days_next, n_days, suffix);
save(fullfile('results', saveStr), ...
    'training_result', 'validation_result', ...
    'Xtrain_norm', 'Ytrain_norm', ...
    'normalized_fields', 'normparams', 'data_test');

% random
training_result = [];
training_result.stepsahead = stepsahead;
training_result.model_inputs = model_inputs;
training_result.model_target = model_target;
training_result.model_excepts = model_excepts;
training_result.hyp = model_random.hyp;
training_result.model = model_random;
Xtrain_norm = X_chosen_random;
Ytrain_norm = y_chosen_random;

validation_result = [];
validation_result.loss_ae = aer;
validation_result.loss_se = ser;
validation_result.loss_lpd = lpdr;
validation_result.loss_mrse = mrser;
validation_result.loss_smse = smser;
validation_result.loss_msll = msllr;
validation_result.loss_rmse = sqrt(ser);

saveStr = sprintf('random_%dinit_%dshift_%dupdate_%dday%s.mat',...
            n_days_init, shift, n_days_next, n_days, suffix);
save(fullfile('results', saveStr), ...
    'training_result', 'validation_result', ...
    'Xtrain_norm', 'Ytrain_norm', ...
    'normalized_fields', 'normparams', 'data_test');

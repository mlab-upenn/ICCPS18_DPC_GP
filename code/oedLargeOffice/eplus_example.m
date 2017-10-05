%% load data

n_samples = 1000;

datafname = 'train-unconstrained-LargeOffice';
data_train = load(fullfile('..', 'data', datafname));
data_train = slice_data(data_train, [1, n_samples]);

% Normalize the data (not all fields)
normalized_fields = {'Ambient', 'Humidity', 'TotalLoad', 'ClgSP', 'LgtSP', 'SupplyAirSP', 'ChwSP'};
[data_train_norm, normparams] = normalize_data(data_train, normalized_fields);
y_train_min =  normparams.TotalLoad.min;
y_train_max = normparams.TotalLoad.max;

datafname = 'test-unconstrained-LargeOffice'; %'test-LargeHotel'; 'test-ramped2-LargeHotel'
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
%{
model_inputs = {...
        'TOD', ...
        'DOW', ...
        {'Ambient', 1:-1:0}, ...
        {'Humidity', 2}, ...
        {'TotalLoad', 3:-1:1}, ...
        {'ClgSP', 2:-1:0}, ... 
        {'LgtSP', [3,0]}, ...
        {'SupplyAirSP', 1:-1:0}, ...
        {'ChwSP', 3:-1:0}};
%}    
model_target = 'TotalLoad';
model_excepts = {'TOD', 'DOW'};
stepsahead = 0;

[X_train_norm, y_train_norm] = construct_data(data_train_norm, model_inputs, model_target, stepsahead, model_excepts);
[X_train, y_train] = construct_data(data_train, model_inputs, model_target, stepsahead, model_excepts);
[X_test_norm, y_test_norm] = construct_data(data_test_norm, model_inputs, model_target, stepsahead, model_excepts);
[~, y_test] = construct_data(data_test, model_inputs, model_target, stepsahead, model_excepts);

%% GP training

model = train_gp_truong(X_train_norm, y_train_norm);


%% GP prediction

% prediction on training data
[mu_train, var_train, muf_train, varf_train] = gp(model.hyp, ...
    model.inference_method, model.mean_function, ...
    model.covariance_function, model.likelihood, ...
    X_train_norm, y_train_norm, X_train_norm);
y_mean_train = postNorm(mu_train, y_train_min, y_train_max);
y_var_train = postNormVar(var_train, y_train_min, y_train_max);

% calculate various errors
loss(y_train, y_mean_train, y_var_train);

% plot results
t = [0:length(y_train)-1]';
f1=figure('Name', 'training');
f1 = plotgp(f1, t, y_train, y_mean_train, sqrt(y_var_train));
axis1 = findobj(f1,'Type','axes');
axis1(2).XLim = [0 n_samples];
axis1(1).XLim = [0 n_samples];
axis1(2).YLim = [0 1.5e6];
axis1(1).YLim = [0 2e5];

% prediction on test data
[mu_test, var_test, muf_test, varf_test] = gp(model.hyp, ...
    model.inference_method, model.mean_function, ...
    model.covariance_function, model.likelihood, ...
    X_train_norm, y_train_norm, X_test_norm);
y_mean_test = postNorm(mu_test, y_train_min, y_train_max);
y_var_test = postNormVar(var_test, y_train_min, y_train_max);

% calculate various errors
loss(y_test, y_mean_test, y_var_test);

% plot results
t = [0:length(y_test)-1]';
f1=figure('Name', 'testing');
f1 = plotgp(f1, t, y_test, y_mean_test, sqrt(y_var_test));
axis1 = findobj(f1,'Type','axes');
axis1(2).XLim = [0 n_samples];
axis1(1).XLim = [0 n_samples];
axis1(2).YLim = [0 1.5e6];
axis1(1).YLim = [0 2e5];

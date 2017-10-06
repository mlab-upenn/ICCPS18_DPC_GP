
n_days_init = 21;
n_days_next = 21;
shift = 30; 
n_days = 16;
suffix = '';  % '_prior' or ''

loadStr = sprintf('evolving_%dinit_%dshift_%dupdate_%dday%s.mat',...
            n_days_init, shift, n_days_next, n_days, suffix);
load(fullfile('results', loadStr))        
model_active = training_result.model;
X_chosen_active = Xtrain_norm;
y_chosen_active = Ytrain_norm;

loadStr = sprintf('random_%dinit_%dshift_%dupdate_%dday%s.mat',...
            n_days_init, shift, n_days_next, n_days, suffix);
load(fullfile('results', loadStr))         
model_random = training_result.model;
X_chosen_random = Xtrain_norm;
y_chosen_random = Ytrain_norm;

y_train_min =  normparams.TotalLoad.min;
y_train_max = normparams.TotalLoad.max;

%% print results for different data sets

model_inputs = training_result.model_inputs;
model_target = training_result.model_target;
model_excepts = training_result.model_excepts;
stepsahead = training_result.stepsahead;

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
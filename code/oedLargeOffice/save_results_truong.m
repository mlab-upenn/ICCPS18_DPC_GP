
% specify after how many samples we need to slice
sim_samples = 96*[3, 4, 7, 10, 14, 21];

for ids = 1:numel(sim_samples)
    
    n_samples = sim_samples(ids);
    results(ids).n_samples = n_samples;
    
    % training data
    
    % RANDOM
    datafnamet = 'random_LargeOffice_uniform_2ramped_3input_28day_20171003_1922.mat';
%     datafnamet = 'random_LargeOffice_prbs_2ramped_3input_28day_20171004_2034.mat';
    
    % OED
%     datafnamet = 'doe_LargeOffice_IG_2ramped_3input_21day_20171003_1804.mat';
%     datafnamet = 'doe_LargeOffice_MV_2ramped_3input_28day_20171003_1817.mat';
    
    data_train = load(fullfile('data', datafnamet));
    data_train = slice_data(data_train, [1, n_samples]);
    
    % Normalize the data (not all fields)
    normalized_fields = {'Ambient', 'Humidity', 'TotalLoad', 'ClgSP', 'LgtSP', 'SupplyAirSP', 'ChwSP'};
    [data_train_norm, normparams] = normalize_data(data_train, normalized_fields);
    y_train_min =  normparams.TotalLoad.min;
    y_train_max = normparams.TotalLoad.max;
    
    % test data
    datafname = 'test-unconstrained-LargeOffice'; %'test-LargeOffice'; 'test-ramped2-LargeOffice'
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
    [X_test, y_test] = construct_data(data_test, model_inputs, model_target, stepsahead, model_excepts);

    data.X_train = X_train;
    data.X_train_norm = X_train_norm;
    data.y_train = y_train;
    data.y_train_norm = y_train_norm;
    data.X_test = X_test;
    data.X_test_norm = X_test_norm;
    data.y_test = y_test;
    data.y_test_norm = y_test_norm;
    data.normparams = normparams;
    results(ids).data = data;
    
    % GP training
    model = train_gp_truong(X_train_norm, y_train_norm);
    results(ids).model = model;

    % GP prediction
    % prediction on training data
    [mu_train, var_train, muf_train, varf_train] = gp(model.hyp, ...
        model.inference_method, model.mean_function, ...
        model.covariance_function, model.likelihood, ...
        X_train_norm, y_train_norm, X_train_norm);
    y_mean_train = postNorm(mu_train, y_train_min, y_train_max);
    y_var_train = postNormVar(var_train, y_train_min, y_train_max);

    % calculate various errors
    [train_results.ae, train_results.se, train_results.lpd, ...
        train_results.mrse, train_results.smse, train_results.msll] = loss(y_train, y_mean_train, y_var_train);
    train_results.rmse = sqrt(train_results.se);
    results(ids).train_results = train_results;

    % prediction on test data
    [mu_test, var_test, muf_test, varf_test] = gp(model.hyp, ...
        model.inference_method, model.mean_function, ...
        model.covariance_function, model.likelihood, ...
        X_train_norm, y_train_norm, X_test_norm);
    y_mean_test = postNorm(mu_test, y_train_min, y_train_max);
    y_var_test = postNormVar(var_test, y_train_min, y_train_max);

    % calculate various errors
    [test_results.ae, test_results.se, test_results.lpd, ...
        test_results.mrse, test_results.smse, test_results.msll] = loss(y_test, y_mean_test, y_var_test);
    test_results.rmse = sqrt(test_results.se);
    results(ids).test_results = test_results;
    
end

saveStr = datafnamet;
save(fullfile('results_truong', saveStr), 'results');

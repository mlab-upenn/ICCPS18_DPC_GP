% Update a GP with new data, offline run

%% Params
orig_gp_file = fullfile('..', 'models', 'LargeOffice', 'random_uniform_ramped2_ahead00_21days_truongisokernel');
new_data_file = fullfile('results', 'random_uniform_ramped2_ahead00_21days_truongisokernel_battery_20171006_0216');
new_gp_file = [orig_gp_file '_evolved'];

num_points = 96*21;
retrain_maxskip = 10;

%% Load the original GP and new data, prepare
ORIGGP = load(orig_gp_file);
NEWDATA = load(new_data_file, ...
    'Ambient', 'ChwSP', 'ClgSP', 'Humidity',  'SupplyAirSP', 'TotalLoad', ...
    'HVACLoad', 'LgtSP', 'Time');
NEWDATA.TOD = NEWDATA.Time(:, 3:4)*[1; 1/60];
NEWDATA.DOW = NEWDATA.Time(:, 6);

% Remove non-vector fields from NEWDATA
flds = fieldnames(NEWDATA);
NEWDATA = rmfield(NEWDATA, flds(cellfun(@(f) ~isvector(NEWDATA.(f)), flds)));

% Normalization
normalized_fields = {'Ambient', 'ChwSP', 'ClgSP', 'Humidity',  'SupplyAirSP', 'TotalLoad'};
NEWDATA = normalize_data(NEWDATA, normalized_fields, ORIGGP.normparams);

% Construct the new inputs and targets to GP
[newinputs, newtargets] = construct_data(NEWDATA, ...
    ORIGGP.training_result.model_inputs, ...
    ORIGGP.training_result.model_target, ...
    ORIGGP.training_result.stepsahead, ...
    ORIGGP.training_result.model_excepts);

%% Run update
new_training = online_batch_update(ORIGGP.training_result.model, ...
    struct('x', [ORIGGP.Xtrain_norm; newinputs], 'y', [ORIGGP.Ytrain_norm; newtargets]), ...
    ORIGGP.training_result.hyp, num_points, ...
    'retrain_maxskip', retrain_maxskip);

% Retrain
ORIGGP.Xtrain_norm = new_training.x;
ORIGGP.Ytrain_norm = new_training.y;

[ORIGGP.training_result.hyp, ORIGGP.training_result.flogtheta, ~] = trainGParx(...
    ORIGGP.training_result.hyp, ORIGGP.training_result.model.inference_method,...
    ORIGGP.training_result.model.mean_function, ...
    ORIGGP.training_result.model.covariance_function, ....
    ORIGGP.training_result.model.likelihood, ...
    ORIGGP.Xtrain_norm, ORIGGP.Ytrain_norm, ...
    @minimize_minfunc, -100);

ORIGGP = rmfield(ORIGGP, 'validation_result');

%% Save
save(new_gp_file, '-struct', 'ORIGGP');

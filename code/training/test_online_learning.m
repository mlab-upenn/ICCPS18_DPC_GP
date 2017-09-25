
load feature_selection/fullset_ahead00_28days_sekernel.mat
data = load('../data/random_sampling_uniform_3input_15day_20170921_0625.mat');

n_samples = 100;

model = struct();
model.init_hyp = training_result.hyp;
model.n_samples = n_samples;
model.model_inputs = training_result.model_inputs;
model.model_target = training_result.model_target;
model.model_stepsahead = training_result.stepsahead;
model.model_excepts = training_result.model_excepts;

sampled_data = online_batch_update(model, data);

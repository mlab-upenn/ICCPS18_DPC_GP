function model = get_default_model(model)
% Returns a complete model where missing fields are set to default values.

% Set default model's functions
assert(isstruct(model));
assert(isfield(model, 'covariance_function'), 'Covariance function must be given.');
if ~isfield(model, 'inference_method')
    model.inference_method = 'infExact';
end
if ~isfield(model, 'mean_function')
    model.mean_function = 'meanZero';
end
if ~isfield(model, 'likelihood')
    model.likelihood = 'likGauss';
end

end


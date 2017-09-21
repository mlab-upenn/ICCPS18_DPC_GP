function valResults = control_validation_gpml(model, hyp, Xtrain, Ytrain, Xtest, Ytest, norm_y_min, norm_y_max)
% Validates a GP model.
% INPUTS:
%   model : the model structure, see CONTROL_TRAIN_GPML.
%   hyp : the hyperparameter structure (after training, e.g., with
%           CONTROL_TRAIN_GPML).
%   Xtrain, Ytrain : original inputs and targets for training.
%   Xtest, Ytest : inputs and targets for validation.
%   norm_y_min, norm_y_max : optional values to de-normalize the output.
%           If they are provided, the predictive outputs will be
%           de-normalized before loss metrics are computed against Ytest.
%
% Regarding normalization, the following must hold true:
%   - Xtrain and Ytrain are normalized if the GP model requires so; in
%       other words, they must be exactly the data used to train the GP
%       model and to perform prediction.
%   - Ytest is always unnormalized.
%   - Xtest is normalized if the GP model requires normalized inputs.
%
% OUTPUTS: a structure of fields
%   .y : the predictive outputs (unnormalized)
%   .s2 : the predictive variance (unnormalized)
%   .y_norm, .s2_norm : in the case the GP's output is normalized
%       (norm_y_min and norm_y_max are provided), these are the normalized
%       / original outputs from the GP model.
%   .loss_rmsq, .loss_ae, .loss_se, .loss_lpd, .loss_mrse, .loss_smse,
%   .loss_msll are some loss measures.

narginchk(4, inf);

nD = size(Xtrain, 2);    % input dimension
assert(size(Xtrain,1) == size(Ytrain,1), 'Xtrain and Ytrain must have the same number of rows.');
assert(size(Xtest, 2) == nD, 'Xtest must have the same dimension as Xtrain.');
assert(size(Xtest,1) == size(Ytest,1), 'Xtest and Ytest must have the same number of rows.');

% Set default model's functions
model = get_default_model(model);

% Validate the hyperparameters
validate_gp(hyp, model.inference_method, model.mean_function, model.covariance_function, model.likelihood, nD);

% Normalization
has_normalization = false;
if nargin > 6
    assert(nargin >= 8, 'Not enough arguments for normalization parameters.');
    assert(isscalar(norm_y_min) && isscalar(norm_y_max) && norm_y_min < norm_y_max, 'Invalid normalization parameters.');
    has_normalization = true;
end

valResults = struct();
[valResults.y, valResults.s2] = gp(hyp, model.inference_method, ...
    model.mean_function, model.covariance_function, model.likelihood, ...
    Xtrain, Ytrain, Xtest);

% De-normalize the data if necessary
if has_normalization
    valResults.y_norm = valResults.y;
    valResults.s2_norm = valResults.s2;
    valResults.y = postNorm(valResults.y_norm, norm_y_min, norm_y_max);
    valResults.s2 = postNormVar(valResults.s2_norm, norm_y_min, norm_y_max);
end

% Compute some loss measures
[valResults.loss_ae, valResults.loss_se, valResults.loss_lpd, valResults.loss_mrse, valResults.loss_smse, valResults.loss_msll] = ...
    loss(Ytest, valResults.y, valResults.s2);
valResults.loss_rmse = sqrt(mean((Ytest - valResults.y).^2));

end

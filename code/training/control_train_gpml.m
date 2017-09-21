function [hyp, flogtheta, model] = control_train_gpml(X, Y, model, hyp0, solver, varargin)
% Train an ARX GP model for control.
% X, Y are training data, obtained from construct_data.m
% model is a structure specifying the GP model for training. Fields:
%          inference_method: a GPML inference method
%                            (optional, default: @infExact)
%             mean_function: a GPML mean function
%                            (optional, default: @meanZero)
%       covariance_function: a GPML covariance function
%                likelihood: a GPML likelihood
%                            (optional, default: @likGauss)
%                     prior: optional GPML prior structure.
% hyp0 is the initial hyperparameter structure.
% solver is the requested optimization solver (optional, default:
% minimize_minfunc).
% varargin are optional solver arguments (default: -100 if minimize_minfunc is the solver).
%
% OUTPUT:
%   hyp : the hyperparameters after training.
%   flogtheta : minus log likelihood for the different runs
%   model : the complete model, where missing fields are set to defaults.

narginchk(3,inf);

nD = size(X, 2);    % input dimension
assert(size(X,1) == size(Y,1), 'X and Y must have the same number of rows.');

% Set default model's functions
model = get_default_model(model);

% Set the inference method based on whether priors are given
if isfield(model, 'prior')
    if ischar(model.inference_method)
        origInf = str2func(model.inference_method);
    else
        origInf = model.inference_method;
    end
    theInfMethod = {@infPrior, origInf, model.prior};
else
    theInfMethod = model.inference_method;
end

% Construct default hyperparameters if not given
if ~exist('hyp0', 'var'), hyp0 = struct(); end
if ~isfield(hyp0, 'cov') || isempty(hyp0.cov)
    hyp0.cov = zeros(gpml_numhyps(model.covariance_function, nD), 1);
end
if ~isfield(hyp0, 'mean') || isempty(hyp0.mean)
    hyp0.mean = zeros(gpml_numhyps(model.mean_function, nD), 1);
end
if ~isfield(hyp0, 'lik') || isempty(hyp0.lik)
    hyp0.lik = zeros(gpml_numhyps(model.likelihood, nD), 1);
end

% Set the solver
if ~exist('solver', 'var') || isempty(solver)
    solver = 'minimize_minfunc';
end
if isempty(varargin)
    if isa(solver, 'function_handle')
        solver = func2str(solver);
    end
    if strcmp(solver, 'minimize_minfunc')
        varargin = {-100};
    end
end
if ischar(solver), solver = str2func(solver); end

% Validate the hyperparameters
validate_gp(hyp0, theInfMethod, model.mean_function, model.covariance_function, model.likelihood, nD);

% Train the model
[hyp, flogtheta, ~] = trainGParx(hyp0, theInfMethod, ...
    model.mean_function, model.covariance_function, model.likelihood, ...
    X, Y, solver, varargin{:});

end
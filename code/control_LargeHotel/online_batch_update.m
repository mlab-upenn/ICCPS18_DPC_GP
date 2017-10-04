function sampled_data = online_batch_update(model, data, init_hyp, n_samples)

% INPUTS:
% data should be struct with fields 'x' and 'y'
%
% model should be struct with fields: covariance_function, mean_function,
% likelihood
%
% init_hyp: initial hyperparams
% n_samples: batch size
% 
% OUTPUT:
% sampled_data is a struct with fields:
%             x: new features for training
%             y: new target values for training

%% Argument checking
assert(all(isfield(model, {'covariance_function', 'mean_function', 'likelihood'})));
assert(isstruct(data) && all(isfield(data, {'x', 'y'})));
assert(n_samples > 0);

%% Validate the hyperparameters
nD = size(data.x, 2);
validate_gp(init_hyp, @infExact, model.mean_function, model.covariance_function, model.likelihood, nD);


%% Convert GPML model to gpml_extension model
model = convert_model(model);


%% Priors on each log covariance parameter
priors.cov = cell(1,numel(init_hyp.cov));
for idc = 1:numel(init_hyp.cov)
    priors.cov{idc}  = get_prior(@gaussian_prior, init_hyp.cov(idc), 1);
%     priors.cov{idc}  = get_prior(@gaussian_prior, 0, 1);
end

% prior on log noise
priors.lik  = {get_prior(@gaussian_prior, init_hyp.lik, 1)};
% priors.lik  = {get_prior(@gaussian_prior, 0, 1)};

% prior on mean hyperparameters
priors.mean = cell(1,numel(init_hyp.mean));
for idc = 1:numel(init_hyp.mean)
    priors.mean{idc}  = get_prior(@gaussian_prior, init_hyp.mean(idc), 1);
%     priors.mean{idc}  = get_prior(@gaussian_prior, 0, 1);
end

model.prior = get_prior(@independent_prior, priors);
model.inference_method = add_prior_to_inference_method(@exact_inference, model.prior);

% setup problem struct
problem = struct;
problem.num_evaluations  = n_samples;
problem.candidate_x_star = data.x;

% function is a simple lookup table
problem.f = @(x,idx) data.y(idx);

% actively learn GP hyperparameters
results = learn_gp_hyperparameters_evolve(problem, model);

sampled_data = struct('x', results.chosen_x, 'y', results.chosen_y);

end


function newmodel = convert_model(model)
% Convert a model from GPML to gpml_extension that can be used for OED
cov_from = {'covSEard', 'covSEiso', 'covMask', ...
    'covProd', 'covSum', ...
    'covConst', 'covRQiso'};
cov_to = {'ard_sqdexp_covariance', 'isotropic_sqdexp_covariance', 'mask_covariance', ...
    'covariance_product', 'covariance_sum', ...
    'constant_covariance', 'iso_rq_covariance'};
mean_from = {'meanConst', 'meanZero'};
mean_to = {'constant_mean', 'zero_mean'};

newmodel = struct;
if isfield(model, 'covariance_function')
    newmodel.covariance_function = iter_convert(model.covariance_function, cov_from, cov_to);
end
    
if isfield(model, 'mean_function')
    newmodel.mean_function = iter_convert(model.mean_function, mean_from, mean_to);
end

if isfield(model, 'likelihood')
    newmodel.likelihood = model.likelihood;
end
end

function newcell = iter_convert(oldcell, from, to)
if iscell(oldcell)
    newcell = cellfun(@(x) iter_convert(x, from, to), oldcell, 'UniformOutput', false);
else
    if isa(oldcell, 'function_handle')
        oldcell = func2str(oldcell);
    end
    
    if ischar(oldcell)
        % Convert function
        [isin, idx] = ismember(oldcell, from);
        if isin
            newcell = str2func(to{idx});
        else
            error('%s is not supported.', oldcell);
        end
    else
        % Keep the original value
        newcell = oldcell;
    end
end
end
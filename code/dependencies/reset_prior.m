function model = reset_prior(model, results)

D = size(results.chosen_x,2); % input space dimension

% covariance function
hyp0.cov = [...
    zeros(1, D), 0, ...
    ]';
covariance_function = {'covSEard'};

% gaussian likelihood function
hyp0.lik = log(0.0005);
likelihood = @likGauss;

% inference method
inference_method = @infExact;

% choose mean function
hyp0.mean = 0;
mean_function = @meanConst;

% solver
solver = @minimize_minfunc;
options = struct('Display', 'off', 'MaxFunEvals', 100);

[new_hyp, ~, ~] = trainGParx(hyp0, inference_method,...
                                    mean_function, covariance_function, ...
                                    likelihood, results.chosen_x, results.chosen_y,...
                                    solver, options);
                                    
% priors on each log covariance parameter
priors.cov = cell(1,D+1);
for idc = 1:D+1
    priors.cov{idc}  = get_prior(@gaussian_prior, new_hyp.cov(idc), 1);
end

% prior on log noise
priors.lik  = {get_prior(@gaussian_prior, new_hyp.lik, 1)};

% prior on constant mean
priors.mean = {get_prior(@gaussian_prior, new_hyp.mean, 1)};

model.prior = get_prior(@independent_prior, priors);
model.inference_method = add_prior_to_inference_method(@exact_inference, model.prior);
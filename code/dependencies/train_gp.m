function model = train_gp(X, y)

D = size(X,2); % input space dimension

% covariance function
hyp0.cov = [...
    zeros(1, D), 0, ...
    ]';
covariance_function = {'covSEard'};
model.covariance_function = covariance_function;

% gaussian likelihood function
hyp0.lik = log(0.0005);
likelihood = @likGauss;
model.likelihood = likelihood;

% inference method
inference_method = @infExact;
model.inference_method = inference_method;

% choose mean function
hyp0.mean = 0;
mean_function = @meanConst;
model.mean_function = mean_function;

% solver
solver = @minimize_minfunc;
options = struct('Display', 'off', 'MaxFunEvals', 100);

[hyp, ~, ~] = trainGParx(hyp0, inference_method,...
    mean_function, covariance_function, ...
    likelihood, X, y,...
    solver, options);

model.hyp0 = hyp0;
model.hyp = hyp;

function model = train_gp_final(X, y)

D = size(X,2); % input space dimension

% covariance function
covariance_function = {'covProd', {{'covSum', {'covConst', {'covMask', {2:D, 'covSEard'}}}}, {'covMask', {1, 'covRQiso'}}}};
hyp0.cov = zeros(1,eval(feval(covariance_function{:})))';
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

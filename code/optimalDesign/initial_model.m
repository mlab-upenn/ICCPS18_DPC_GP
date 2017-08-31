function hyp = initial_model(file, n_samples)

%% load data

ctrl_horizon = 1;
order_autoreg = 3;
plot_results = 0;

[X, y] = load_data(file, order_autoreg, ctrl_horizon);

X_train = X(1:n_samples,:);
y_train = y(1:n_samples,:);
X_test = X(n_samples+1:2*n_samples,:);
y_test = y(n_samples+1:2*n_samples,:);

% standardize the data set
[X_train_norm, X_train_min, X_train_max] = preNorm(X_train);
[y_train_norm, y_train_min, y_train_max] = preNorm(y_train);

X_test_norm = preNorm(X_test, X_train_min, X_train_max);
y_test_norm = preNorm(y_test, y_train_min, y_train_max);

%% GP model definition

D = size(X,2); % input space dimension

% covariance function
cov = {'covSEard'};
hyp0.cov = [...
    zeros(1, D), 0, ...
    ]';

% gaussian likelihood function
lik = @likGauss;
hyp0.lik = log(0.01);

% inference method
inf = @infExact;

% choose mean function
meanf = @meanConst;
hyp0.mean = 0;

% solver
solver = @minimize_minfunc;
maxiter = -100;

% training
[hyp, flogtheta, ~] = trainGParx(hyp0, inf, meanf, cov, lik, X_train_norm, y_train_norm, solver, -100);

% prediction on training data
[mu_train, var_train, muf_train, varf_train] = gp(hyp, inf, meanf, cov, lik, X_train_norm, y_train_norm, X_train_norm);
y_mean_train = postNorm(mu_train, y_train_min, y_train_max);
y_var_train = postNormVar(var_train, y_train_min, y_train_max);

% calculate various errors
loss(y_train, y_mean_train, y_var_train);

% plot results
if plot_results
    t = [0:length(y_train)-1]';
    f1 = figure;
    title('training');
    f1 = plotgp(f1, t, y_train, y_mean_train, sqrt(y_var_train));
    axis1 = findobj(f1,'Type','axes');
    axis1(2).XLim = [0 n_samples];
    axis1(1).XLim = [0 n_samples];
    axis1(1).YLim = [0 5e4];
    axis1(2).YLim = [0 5e5];
end

% prediction on test data
[mu_test, var_test, muf_test, varf_test] = gp(hyp, inf, meanf, cov, lik, X_train_norm, y_train_norm, X_test_norm);
y_mean_test = postNorm(mu_test, y_train_min, y_train_max);
y_var_test = postNormVar(var_test, y_train_min, y_train_max);

% calculate various errors
loss(y_test, y_mean_test, y_var_test);

% plot results
if plot_results
    t = [0:length(y_test)-1]';
    f2 = figure;
    title('testing');
    f2 = plotgp(f2, t, y_test, y_mean_test, sqrt(y_var_test));
    axis1 = findobj(f2,'Type','axes');
    axis1(2).XLim = [0 n_samples];
    axis1(1).XLim = [0 n_samples];
    axis1(1).YLim = [0 5e4];
    axis1(2).YLim = [0 5e5];
end

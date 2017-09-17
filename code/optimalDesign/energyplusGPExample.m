%% load data

file = 'unconstrained-LargeHotel';
ctrl_horizon = 1;
order_autoreg = 3;
n_samples = 1000;

% ctrl_variables = {'ClgSP', 'KitchenClgSP', 'GuestClgSP', 'SupplyAirSP', 'ChwSP'};
ctrl_vars = {'GuestClgSP', 'SupplyAirSP', 'ChwSP'};

[X, y] = load_data(file, order_autoreg, ctrl_vars);
% standardize the data set
[X_train_norm, X_train_min, X_train_max] = preNorm(X);
[y_train_norm, y_train_min, y_train_max] = preNorm(y);

X_train = X(1:n_samples,:);
y_train = y(1:n_samples);
X_train_norm = X_train_norm(1:n_samples,:);
y_train_norm = y_train_norm(1:n_samples,:);

% X_test = X(n_samples+1:2*n_samples,:);
% y_test = y(n_samples+1:2*n_samples,:);
datafile = 'test-LargeHotel';
[X_test, y_test] = load_data(datafile, order_autoreg, ctrl_vars);
X_test_norm = preNorm(X_test, X_train_min, X_train_max);
y_test_norm = preNorm(y_test, y_train_min, y_train_max);

%% GP model definition

D = size(X,2); % input space dimension

% covariance function
cov = {'covSEard'};
hyp0.cov = [...
    zeros(1, D), 0, ...
    ]';
% cov = {'covProd', {...
%     {'covSum', {'covConst', 'covSEard'}},...
%     {'covRQiso'}} ...  % temporal
%     };
% hyp0.cov = [...
%     0, ...
%     zeros(1, D), 0, ...    % non-temporal
%     0, 0, 0, ...  % temporal
%     ]';

% gaussian likelihood function
lik = @likGauss;
% prior.lik = {{@priorDelta}};
hyp0.lik = log(0.0005);

% inference method
inf = @infExact;
% inf = {@infPrior,@infExact,prior};

% choose mean function
meanf = @meanConst;
hyp0.mean = 0;

% solver
solver = @minimize_minfunc;
maxiter = -100;

% training
[hyp, flogtheta, ~] = trainGParx(hyp0, inf, meanf, cov, lik, X_train_norm, y_train_norm, solver, maxiter);

% prediction on training data
[mu_train, var_train, muf_train, varf_train] = gp(hyp, inf, meanf, cov, lik, X_train_norm, y_train_norm, X_train_norm);
y_mean_train = postNorm(mu_train, y_train_min, y_train_max);
y_var_train = postNormVar(var_train, y_train_min, y_train_max);

% calculate various errors
loss(y_train, y_mean_train, y_var_train);

% plot results
t = [0:length(y_train)-1]';
f1=figure('Name', 'training');
f1 = plotgp(f1, t, y_train, y_mean_train, sqrt(y_var_train));
axis1 = findobj(f1,'Type','axes');
axis1(2).XLim = [0 n_samples];
axis1(1).XLim = [0 n_samples];
axis1(1).YLim = [0 5e4];
axis1(2).YLim = [0 5e5];

% prediction on test data
[mu_test, var_test, muf_test, varf_test] = gp(hyp, inf, meanf, cov, lik, X_train_norm, y_train_norm, X_test_norm);
y_mean_test = postNorm(mu_test, y_train_min, y_train_max);
y_var_test = postNormVar(var_test, y_train_min, y_train_max);

% calculate various errors
loss(y_test, y_mean_test, y_var_test);

% plot results
t = [0:length(y_test)-1]';
f1=figure('Name', 'testing');
f1 = plotgp(f1, t, y_test, y_mean_test, sqrt(y_var_test));
axis1 = findobj(f1,'Type','axes');
axis1(2).XLim = [0 n_samples];
axis1(1).XLim = [0 n_samples];
axis1(1).YLim = [0 5e4];
axis1(2).YLim = [0 5e5];

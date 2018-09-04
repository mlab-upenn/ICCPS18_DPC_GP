%% predict on non-optimial hyperparams

rng(1);

% Test data
n = 100;
Xtest = linspace(-5, 5, n)';

param = 0.1;
K_ss = kernel(Xtest, Xtest, param);

% Noiseless training data
Xtrain = [-4, -3, -2, -1, 1]';
ytrain = [-2, 0, 1, 2, -1]';

% Apply the kernel function to our training points
K = kernel(Xtrain, Xtrain, param);
L = chol(K + 0.00005*eye(numel(Xtrain)));

% Compute the mean at our test points.
K_s = kernel(Xtrain, Xtest, param);
Lk = L\K_s;
mu = Lk' * (L\ytrain);

% Compute the standard deviation so we can plot it
s2 = diag(K_ss) - sum(Lk.^2, 1)';
stdv = sqrt(s2);

% Draw samples from the posterior at our test points.
L = chol(K_ss + 1e-6*eye(n) - (Lk'*Lk));
f_post = repmat(mu, [1,3]) + L * mvnrnd(zeros(n,1), eye(n), 3)';

D = size(Xtest,2);
hyp0.cov = [...
    log(sqrt(param)), 0, ...
    ]';
cov = {'covSEiso'};
hyp0.lik = log(sqrt(0.00005));
lik = @likGauss;
inf = @infExact;
mean = @meanZero;
[mugp, s2gp] = gp(hyp0, inf, mean, cov, lik, Xtrain, ytrain, Xtest);
stdvgp = sqrt(s2gp);

figure(); hold on;
xfill = [Xtest; flipdim(Xtest,1)]; 
yfill = [mugp+2*stdvgp; flipdim(mugp-2*stdvgp,1)];
fill(xfill, yfill, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
alpha(0.75);
plot(Xtest, f_post(:,1), 'r', 'LineWidth', 1);
plot(Xtest, f_post(:,2), 'b', 'LineWidth', 1);
plot(Xtest, f_post(:,3), 'g', 'LineWidth', 1);
plot(Xtest, mugp, 'k--', 'LineWidth', 2);
plot(Xtrain, ytrain, 'kp', 'MarkerSize', 15, 'MarkerFaceColor', 'k');
axis([-5, 5, -3, 3])
title('Three samples from the GP posterior')
% http://katbailey.github.io/post/gaussian-processes-for-dummies/

rng(1);

% Test data
n = 100;
Xtest = linspace(-5, 5, n)';

param = 0.1;
K_ss = kernel(Xtest, Xtest, param);
mu = zeros(n,1);
s2 = diag(K_ss);
stdv = sqrt(s2);

% Get cholesky decomposition (square root) of the
% covariance matrix
L = chol(K_ss + 1e-10*eye(n));

% Sample 3 sets of standard normals for our test points,
% multiply them by the square root of the covariance matrix
f_prior = L * mvnrnd(zeros(n,1), eye(n), 3)';

% Now let's plot the 3 sampled functions.
figure(); hold on; reset_color;
plot(Xtest, mu, 'k--', 'LineWidth', 2);
xfill = [Xtest; flipdim(Xtest,1)]; 
yfill = [mu+2*stdv; flipdim(mu-2*stdv,1)];
fill(xfill, yfill, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
alpha(0.75);
plot(Xtest, f_prior(:,1), 'r', 'LineWidth', 1);
plot(Xtest, f_prior(:,2), 'b', 'LineWidth', 1);
plot(Xtest, f_prior(:,3), 'g', 'LineWidth', 1);
axis([-5, 5, -3, 3])
title('Three samples from the GP prior')

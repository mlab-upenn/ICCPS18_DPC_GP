%% predict on optimial hyperparams

plot_settings;
rng(1);

% Test data
n = 100;
Xtest = linspace(-4, 4, n)';

param = 0.1;

% Noiseless training data
% Xtrain = [-3, -2, -1, 0, 2]';
% ytrain = [-2, 0, 1, 2, -1]';

Xtrain = [-3, -1, -1.5, 4]';
myfunc = @(x) sin(x);
ytrain = myfunc(Xtrain);

D = size(Xtest,2);
hyp0.cov = [...
    log(sqrt(param)), 0, ...
    ]';
cov = {'covSEiso'};
hyp0.lik = log(sqrt(0.00005));
lik = @likGauss;
inf = @infExact;
mean = @meanZero;

solver = @minimize_minfunc;
options = struct('Display', 'on', 'MaxFunEvals', 100);
[hyp, ~, ~] = trainGParx(hyp0, inf, mean, cov, lik, Xtrain, ytrain,...
    solver, options);
[mugp, s2gp] = gp(hyp, inf, mean, cov, lik, Xtrain, ytrain, Xtest);
stdvgp = sqrt(s2gp);

figure; hold on;
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 4];
xfill = [Xtest; flipdim(Xtest,1)]; 
yfill = [mugp+2*stdvgp; flipdim(mugp-2*stdvgp,1)];
h1 = fill(xfill, yfill, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
alpha(0.75);
h2 = plot(Xtest, mugp, 'k--', 'LineWidth', 2);
h3 = plot(Xtrain, ytrain, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
legend([h1, h2, h3], '\mu \pm 2\sigma', '\mu', 'training sample', 'Location', 'North', 'Orientation', 'horizontal');
axis([-4, 4, -3, 3])
xlabel('input x')
ylabel('output y')
print('prior', '-dpng')

figure; hold on;
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2];
xfill = [Xtest; flipdim(Xtest,1)]; 
yfill = [stdvgp.^2; 0*flipdim(stdvgp,1)];
fill(xfill, yfill, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
alpha(0.75);
axis([-4, 4, 0, 0.5])
xlabel('input x')
ylabel('utility')
print('prior-util', '-dpng')

% update based on maximum variance
[~, ind] = max(s2gp);
Xtrain = [Xtrain; Xtest(ind)];
ytrain = myfunc(Xtrain);

[hyp, ~, ~] = trainGParx(hyp0, inf, mean, cov, lik, Xtrain, ytrain,...
    solver, options);
[mugp, s2gp] = gp(hyp, inf, mean, cov, lik, Xtrain, ytrain, Xtest);
stdvgp = sqrt(s2gp);

figure; hold on;
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 4];
xfill = [Xtest; flipdim(Xtest,1)]; 
yfill = [mugp+2*stdvgp; flipdim(mugp-2*stdvgp,1)];
h1 = fill(xfill, yfill, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
alpha(0.75);
h2 = plot(Xtest, mugp, 'k--', 'LineWidth', 2);
h3 = plot(Xtrain, ytrain, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r');
legend([h1, h2, h3], '\mu \pm 2\sigma', '\mu', 'training sample', 'Location', 'North', 'Orientation', 'horizontal');
axis([-4, 4, -3, 3])
xlabel('input x')
ylabel('output y')
print('posterior', '-dpng', '-r600')

figure; hold on;
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 2];
xfill = [Xtest; flipdim(Xtest,1)]; 
yfill = [stdvgp.^2; 0*flipdim(stdvgp,1)];
fill(xfill, yfill, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
alpha(0.75);
axis([-4, 4, 0, 0.5])
xlabel('input x')
ylabel('utility')
print('posterior-util', '-dpng', '-r600')

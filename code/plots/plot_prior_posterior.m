
load ../training/feature_selection_LargeHotel/doe_sampling_noreset_IG_2ramped_3input_ahead00_14days_truongkernel.mat

model = training_result.model;
hyp = training_result.hyp;
cov = model.covariance_function;
lik = model.likelihood;
inff = model.inference_method;
meanf = model.mean_function;

n_test_samples = 0.5*96;

y_train_min = normparams.TotalLoad.min;
y_train_max = normparams.TotalLoad.max;
X_sample = Xtrain_norm(1:n_test_samples,:);

% d_n = exp(2*hyp.lik);
% m_ = hyp.mean*ones(1,n_samples);
% m_ss = hyp.mean*ones(1,n_test_samples);
K_ss = feval(cov{:},hyp.cov, X_sample);
% K_s = feval(cov{:},hyp.cov, X_sample, X_train_norm);
% K = feval(cov{:},hyp.cov, X_train_norm);
% 
% n_trajectories = 5;
% y_traj = mvnrnd(m_ss, K_ss, n_trajectories);
% for idt = 1:n_trajectories 
%     y_traj(idt,:) = postNorm(y_traj(idt,:), y_train_min, y_train_max)/1e3;
% end

% prediction with priors
t = [1:n_test_samples]';
y = postNorm(zeros(n_test_samples,1), y_train_min, y_train_max)/1e3-70;
pm = y(1);

std = sqrt(postNormVar(diag(K_ss), y_train_min, y_train_max))/1e3/3;
ps = std(1);

ix_plot = 1:length(t);
xfill = [t(ix_plot); flipdim(t(ix_plot),1)]; 
yfill = [y(ix_plot)+2*std(ix_plot);flipdim(y(ix_plot)-2*std(ix_plot),1)]; 

figure; hold on; grid off;
h1 = fill(xfill, yfill, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
h2 = plot(t,y, '--k', 'LineWidth',1); 

% prediction with posterior
[mu, var, muf, varf] = gp(hyp, inff, meanf, cov, lik, Xtrain_norm, Ytrain_norm, X_sample);
y_mean = postNorm(mu, y_train_min, y_train_max);
y_var = postNormVar(var, y_train_min, y_train_max);

t = [1:n_test_samples]';
y = y_mean/1e3;
std = sqrt(y_var)/1e3;

ix_plot = 1:length(t);
xfill = [t(ix_plot); flipdim(t(ix_plot),1)]; 
yfill = [y(ix_plot)+2*std(ix_plot);flipdim(y(ix_plot)-2*std(ix_plot),1)]; 
h3 = fill(xfill, yfill, [249, 229, 255]/255, 'EdgeColor', [249, 229, 255]/255);
h4 = plot(t,y, '-r', 'LineWidth', 0.5); 
ax = gca;
ax.XTick = 0.5*[16, 32, 48, 64, 80];
ax.XTickLabel = {'2am','4am','6am','8am', '10am'}';
axis([1, n_test_samples, pm-2.2*ps,pm+2.2*ps]);
hleg = legend([h2, h1, h4, h3], ['prior \mu'], ['prior \mu \pm 2\sigma'], 'posterior \mu', 'posterior \mu \pm 2\sigma');
set(hleg, 'Location','NorthOutside', 'box', 'off', 'Orientation', 'Horizontal'); 

ylabel('power [kW]')
box on;

cleanfigure()

matlab2tikz('width', '\fwidth', 'height', '\hwidth', '../../paper/src/figures/gp-prior-post.tex',...
    'extraaxisoptions',['xlabel style={font=\footnotesize},'...
                       'ylabel style={font=\footnotesize},',...
                       'legend style={font=\footnotesize},',...
                       'ticklabel style={font=\footnotesize},'...
                       'ylabel shift = -5 pt,'...
                       'xlabel shift = -5 pt,']);

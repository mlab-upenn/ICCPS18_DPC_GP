
load('../optimalDesign/results/selection_20171001_0054.mat');

lw = 0.5;

n_days = 1;
idx = 1:n_days*96;
t = [0:1:length(idx)-1]';

% active learning
figure;
sys = y_test(idx)/1e3;
y = f_star_mean_active(idx)/1e3;
std = sqrt(f_star_variance_active(idx))/1e3;

ix_plot = 1:length(t);
xfill = [t(ix_plot); flipdim(t(ix_plot),1)]; 
yfill = [y(ix_plot)+2*std(ix_plot);flipdim(y(ix_plot)-2*std(ix_plot),1)]; 

h1 = subplot(3,1,1:2); hold on; grid on;
fill(xfill, yfill, [249, 229, 255]/255, 'EdgeColor', [249, 229, 255]/255);
plot(t,y, '--k', 'LineWidth', lw); 
h = plot(t,sys, '-', 'LineWidth', lw);
ax = gca;
ax.XTick = [16, 32, 48, 64, 80];
ax.XTickLabel = {}';
xlim([0 t(end)]);
ylim([0 450]);
ylabel('power [kW]'); 
hold off;
title('OPTIMAL SELECTION')
hleg = legend('\mu \pm 2\sigma', '\mu', 'system','Location','NorthOutside', 'Orientation', 'Horizontal'); 
set(hleg, 'box', 'off');
box on;
 
h2 = subplot(3,1,3); hold on; grid on;
xfill = [t(ix_plot); flipdim(t(ix_plot),1)]; 
yfill = [2*std(ix_plot);zeros(size(std(ix_plot)))]; 
fill(xfill, yfill, [249, 229, 255]/255, 'EdgeColor', [249, 229, 255]/255);
plot(t,abs(y-sys), '-', 'LineWidth', lw);
ax = gca;
ax.YTick = [0, 40, 80];
ax.XTick = [16, 32, 48, 64, 80];
ax.XTickLabel = {'4am','8am','12pm','4pm', '8pm'}';
xlim([0 t(end)]);
ylabel('error [kW]'); 
% xlabel('time');
box on;

cleanfigure()
matlab2tikz('width', '\fwidth', 'height', '\hwidth', '../../paper/src/figures/batch-oed.tex',...
    'extraaxisoptions',['xlabel style={font=\footnotesize},'...
                       'ylabel style={font=\footnotesize},',...
                       'legend style={font=\footnotesize},',...
                       'ticklabel style={font=\footnotesize},'...
                       'title style={font=\normalsize},'...
                       'title style={yshift=2ex},'...                       
                       'ylabel shift = -5 pt,'...
                       'xlabel shift = -5 pt,']);

% random sampling
figure;
sys = y_test(idx)/1e3;
y = f_star_mean(idx)/1e3;
std = sqrt(f_star_variance(idx))/1e3;

ix_plot = 1:length(t);
xfill = [t(ix_plot); flipdim(t(ix_plot),1)]; 
yfill = [y(ix_plot)+2*std(ix_plot);flipdim(y(ix_plot)-2*std(ix_plot),1)]; 

h1 = subplot(3,1,1:2); hold on; grid on;
fill(xfill, yfill, [249, 229, 255]/255, 'EdgeColor', [249, 229, 255]/255);
plot(t,y, '--k', 'LineWidth', lw); 
plot(t,sys, '-', 'LineWidth', lw);
xlim([0 t(end)]);
ylim([0 450]);
ax = gca;
ax.XTick = [16, 32, 48, 64, 80];
ax.XTickLabel = {}';
ylabel('power [kW]');
hold off;
title('RANDOM SAMPLING')
hleg = legend('\mu \pm 2\sigma', '\mu', 'system','Location','NorthOutside', 'Orientation', 'Horizontal'); 
set(hleg, 'box', 'off');
box on;
 
h2 = subplot(3,1,3); hold on; grid on;
xfill = [t(ix_plot); flipdim(t(ix_plot),1)]; 
yfill = [2*std(ix_plot);zeros(size(std(ix_plot)))]; 
fill(xfill, yfill, [249, 229, 255]/255, 'EdgeColor', [249, 229, 255]/255);
plot(t,abs(y-sys), '-', 'LineWidth', lw);
ax = gca;
ax.YTick = [0, 40, 80];
ax.XTick = [16, 32, 48, 64, 80];
ax.XTickLabel = {'4am','8am','12pm','4pm', '8pm'}';
xlim([0 t(end)]);
ylabel('error [kW]'); 
% xlabel('time');
box on;

cleanfigure()
matlab2tikz('width', '\fwidth', 'height', '\hwidth', '../../paper/src/figures/batch-rand.tex',...
    'extraaxisoptions',['xlabel style={font=\footnotesize},'...
                       'ylabel style={font=\footnotesize},',...
                       'legend style={font=\footnotesize},',...
                       'ticklabel style={font=\footnotesize},'...
                       'title style={font=\normalsize},'...
                       'title style={yshift=2ex},'...
                       'ylabel shift = -5 pt,'...
                       'xlabel shift = -5 pt,']);

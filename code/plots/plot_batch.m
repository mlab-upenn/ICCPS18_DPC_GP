
load('../optimalDesign/results/selection_20171001_0054.mat');

n_days = 1;
idx = 1:n_days*96;
t = [0:1:length(idx)-1]';

figure;

% active learning
sys = y_test(idx)/1e3;
y = f_star_mean_active(idx)/1e3;
std = sqrt(f_star_variance_active(idx))/1e3;

ix_plot = 1:length(t);
xfill = [t(ix_plot); flipdim(t(ix_plot),1)]; 
yfill = [y(ix_plot)+2*std(ix_plot);flipdim(y(ix_plot)-2*std(ix_plot),1)]; 

h1 = subplot(2,1,1); hold on; grid on; box on;
fill(xfill, yfill, [249, 229, 255]/255, 'EdgeColor', [249, 229, 255]/255);
plot(t,y, '-', 'LineWidth',1); 
h = plot(t,sys, '-', 'LineWidth',1);
ax = gca;
ax.XTick = [16, 32, 48, 64, 80];
ax.XTickLabel = {'4am','8am','12pm','4pm', '8pm'}';
xlim([0 t(end)]);
ylim([0 450]);
ylabel('power [kW]'); 
hleg = legend('\mu \pm 2\sigma', '\mu', 'system','Location','NorthOutside', 'Orientation', 'Horizontal'); 
set(hleg, 'box', 'off');

% random sampling
sys = y_test(idx)/1e3;
y = f_star_mean(idx)/1e3;
std = sqrt(f_star_variance(idx))/1e3;

ix_plot = 1:length(t);
xfill = [t(ix_plot); flipdim(t(ix_plot),1)]; 
yfill = [y(ix_plot)+2*std(ix_plot);flipdim(y(ix_plot)-2*std(ix_plot),1)]; 

h2 = subplot(2,1,2); hold on; grid on; box on;
h1.Position(4) = h2.Position(4);
fill(xfill, yfill, [249, 229, 255]/255, 'EdgeColor', [249, 229, 255]/255);
plot(t,y, '-', 'LineWidth',1); 
plot(t,sys, '-', 'LineWidth',1);
xlim([0 t(end)]);
ylim([0 450]);
ax = gca;
ax.XTick = [16, 32, 48, 64, 80];
ax.XTickLabel = {'4am','8am','12pm','4pm', '8pm'}';
% xlabel('time'); 
ylabel('power [kW]');  
% legend('\mu \pm 2\sigma', '\mu', 'system','Location','NorthEast');
cleanfigure()

matlab2tikz('width', '\fwidth', 'height', '\hwidth', '../../paper/src/figures/batch-acc.tex',...
    'extraaxisoptions',['xlabel style={font=\footnotesize},'...
                       'ylabel style={font=\footnotesize},',...
                       'legend style={font=\footnotesize},',...
                       'ticklabel style={font=\footnotesize},'...
                       'ylabel shift = -5 pt,'...
                       'xlabel shift = -5 pt,']);

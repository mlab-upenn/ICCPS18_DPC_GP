
load('selection_20170930_0037.mat')

t = [0:length(y_test)-1]';
sys = y_test/1e3;
y = f_star_mean_active/1e3;
std = sqrt(f_star_variance_active)/1e3;

ix_plot = 1:length(t);
xfill = [t(ix_plot); flipdim(t(ix_plot),1)]; 
yfill = [y(ix_plot)+2*std(ix_plot);flipdim(y(ix_plot)-2*std(ix_plot),1)]; 

figure; 
h1 = subplot(2,1,1); hold on; grid on; box on;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 6 3];

fill(xfill, yfill, [249, 229, 255]/255, 'EdgeColor', [249, 229, 255]/255);
plot(t,y, '-', 'LineWidth',1); 
plot(t,sys, '-', 'LineWidth',1);
xlim([0 100]);
ylim([0 450]);
ylabel('power [kW]'); 
hleg = legend('\mu \pm 2\sigma', '\mu', 'system','Location','NorthOutside', 'Orientation', 'Horizontal'); 
set(hleg, 'box', 'off');
% matlab2tikz('width', '\fwidth', 'height', '\hwidth', '../../../paper/src/figures/selection_active.tex');

t = [0:length(y_test)-1]';
sys = y_test/1e3;
y = f_star_mean/1e3;
std = sqrt(f_star_variance)/1e3;

ix_plot = 1:length(t);
xfill = [t(ix_plot); flipdim(t(ix_plot),1)]; 
yfill = [y(ix_plot)+2*std(ix_plot);flipdim(y(ix_plot)-2*std(ix_plot),1)]; 

h2 = subplot(2,1,2); hold on; grid on; box on;
fill(xfill, yfill, [249, 229, 255]/255, 'EdgeColor', [249, 229, 255]/255);
plot(t,y, '-', 'LineWidth',1); 
plot(t,sys, '-', 'LineWidth',1); 
xlim([0 100]);
ylim([0 450]);
xlabel('sample index [-]'); 
ylabel('power [kW]');  
% legend('\mu \pm 2\sigma', '\mu', 'system','Location','NorthEast');
cleanfigure()
% matlab2tikz('width', '\fwidth', 'height', '\hwidth', '../../../paper/src/figures/selection_random.tex');
matlab2tikz('width', '\fwidth', 'height', '\hwidth', '../../../paper/src/figures/selection.tex');





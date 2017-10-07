
% figure settings
lw = 0.75;
ms = 6;

load('random_uniform_ramped2_ahead00_21days_truongisokernel_battery_20171005_0500.mat')

shift = 6;
idx = 96*shift + (1:96);
t = 1:96;

battery = BatteryPower(idx)/1e6;
system = TotalLoad(idx)/1e6;
baseline = baseline(idx)/1e6;
mu = PowerExpected(idx)/1e6;
std = sqrt(PowerVariance(idx))/1e6;
reference = ref(idx)/1e6;
ClgSP = ClgSP(idx);
SupplyAirSP = SupplyAirSP(idx);
ChwSP = ChwSP(idx);

pstart = 7*4;
pend = 20*4;
pidx = pstart:pend;

% power consumption
figure; hold on; grid on; box on;
h1 = plot(t(pidx), system(pidx), '-', 'LineWidth', lw);
h2 = plot(t(pidx), baseline(pidx), '-', 'LineWidth', lw, 'Color', [0.4940    0.1840    0.5560]);
h3 = plot(t(pidx), mu(pidx), '-', 'LineWidth', lw, 'Color', [0.8500    0.3250    0.0980]);
h4 = plot(t(pidx), reference(pidx), '--', 'LineWidth', lw, 'Color', [0.4660    0.6740    0.1880]);

plot(52*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);
plot(56*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);
plot(64*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);
plot(68*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);

axis([t(pidx(1)), t(pidx(end)), 900/1e3, 1400/1e3])
ax = gca;
ax.XTick = [40, 52, 56, 64, 68];
ax.XTickLabel = {'10am', '1pm', '2pm', '4pm', '5pm'}';
ylabel('power [MW]')

hleg = legend([h1, h2, h3, h4], 'system', 'baseline', '\mu', 'reference');
set(hleg, 'box', 'off', 'Location', 'NorthOutside', 'Orientation', 'Horizontal');

matlab2tikz('width', '\fwidth', 'height', '\hwidth', '../../paper/src/figures/control-tracking.tex',...
    'extraaxisoptions',['xlabel style={font=\footnotesize},'...
                       'ylabel style={font=\footnotesize},',...
                       'legend style={font=\footnotesize},',...
                       'ticklabel style={font=\footnotesize},'...
                       'ylabel shift = -5 pt,'...
                       'xlabel shift = -5 pt,']);

% control inputs
figure; hold on; grid on; box on;
h1 = plot(t(pidx), ClgSP(pidx), '-', 'LineWidth', lw);
h2 = plot(t(pidx), SupplyAirSP(pidx), '-', 'LineWidth', lw, 'Color', [0.4940    0.1840    0.5560]);
h3 = plot(t(pidx), ChwSP(pidx), '-', 'LineWidth', lw, 'Color', [0.8500    0.3250    0.0980]);

plot(52*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);
plot(56*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);
plot(64*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);
plot(68*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);

axis([t(pidx(1)), t(pidx(end)), 5, 30])
ax = gca;
ax.XTick = [40, 52, 56, 64, 68];
ax.XTickLabel = {'10am', '1pm', '2pm', '4pm', '5pm'}';
ylabel('temperature [^oC]')

hleg = legend([h1, h2, h3], 'cooling', 'supply air', 'chilled water');
set(hleg, 'box', 'on', 'Location', 'west', 'Orientation', 'Horizontal', 'EdgeColor',[1 1 1]);

matlab2tikz('width', '\fwidth', 'height', '\hwidth', '../../paper/src/figures/control-inputs.tex',...
    'extraaxisoptions',['xlabel style={font=\footnotesize},'...
                       'ylabel style={font=\footnotesize},',...
                       'legend style={font=\footnotesize},',...
                       'ticklabel style={font=\footnotesize},'...
                       'ylabel shift = -5 pt,'...
                       'xlabel shift = -5 pt,']);
                   
% error and std
figure; hold on; box on;
t = t(pidx)';
mu = mu(pidx)*1e3;
std = std(pidx)*1e3;
sys = system(pidx)*1e3;
err = abs(mu-sys);
err(t<52 | t>68) = NaN;

ix_plot = 1:length(t);
xfill = [t(ix_plot); flipdim(t(ix_plot),1)]; 
yfill = [2*std(ix_plot);zeros(size(std(ix_plot)))]; 
h1 = fill(xfill, yfill, [249, 229, 255]/255, 'EdgeColor', [249, 229, 255]/255);
h2 = plot(t, err, '-', 'LineWidth', lw);

grid on;

plot(52*ones(1,10), linspace(0,150,10), '--k', 'LineWidth', lw);
plot(56*ones(1,10), linspace(0,150,10), '--k', 'LineWidth', lw);
plot(64*ones(1,10), linspace(0,150,10), '--k', 'LineWidth', lw);
plot(68*ones(1,10), linspace(0,150,10), '--k', 'LineWidth', lw);

axis([t(1), t(end), 0, 150])
ax = gca;
ax.XTick = [40, 52, 56, 64, 68];
ax.XTickLabel = {'10am', '1pm', '2pm', '4pm', '5pm'}';
ylabel('error [kW]')

hleg = legend([h2, h1], 'prediction error', '2\sigma');
set(hleg, 'box', 'on', 'Location', 'north', 'Orientation', 'Horizontal', 'EdgeColor',[1 1 1]);

matlab2tikz('width', '\fwidth', 'height', '\hwidth', '../../paper/src/figures/control-error.tex',...
    'extraaxisoptions',['xlabel style={font=\footnotesize},'...
                       'ylabel style={font=\footnotesize},',...
                       'legend style={font=\footnotesize},',...
                       'ticklabel style={font=\footnotesize},'...
                       'ylabel shift = -5 pt,'...
                       'xlabel shift = -5 pt,']);
                   
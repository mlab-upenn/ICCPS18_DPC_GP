
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


% control inputs
figure;

subplot(3,1,1); hold on; grid on; box on;
h1 = plot(t(pidx), ClgSP(pidx), '-', 'LineWidth', lw);

plot(52*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);
plot(56*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);
plot(64*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);
plot(68*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);

axis([t(pidx(1)), t(pidx(end)), 23, 27])
ax = gca;
ax.XTick = [40, 52, 56, 64, 68];
ax.XTickLabel = {'10am', '1pm', '2pm', '4pm', '5pm'}';
ax.YTick = [24, 26];
ylabel('temp. [^oC]')

hleg = legend(h1, 'cooling');
set(hleg, 'box', 'on', 'Location', 'west', 'Orientation', 'Horizontal', 'EdgeColor',[1 1 1]);

% matlab2tikz('width', '\fwidth', 'height', '\hwidth', '../../paper/src/figures/control-cooling.tex',...
%     'extraaxisoptions',['xlabel style={font=\footnotesize},'...
%                        'ylabel style={font=\footnotesize},',...
%                        'legend style={font=\footnotesize},',...
%                        'ticklabel style={font=\footnotesize},'...
%                        'ylabel shift = -5 pt,'...
%                        'xlabel shift = -5 pt,']);

% control inputs
subplot(3,1,2); hold on; grid on; box on;
h2 = plot(t(pidx), SupplyAirSP(pidx), '-', 'LineWidth', lw, 'Color', [0.4940    0.1840    0.5560]);

plot(52*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);
plot(56*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);
plot(64*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);
plot(68*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);

axis([t(pidx(1)), t(pidx(end)), 12.4, 14])
ax = gca;
ax.XTick = [40, 52, 56, 64, 68];
ax.XTickLabel = {'10am', '1pm', '2pm', '4pm', '5pm'}';
ax.YTick = [13, 14];
ylabel('temp. [^oC]')

hleg = legend(h2, 'supply air');
set(hleg, 'box', 'on', 'Location', 'northwest', 'Orientation', 'Horizontal', 'EdgeColor',[1 1 1]);

% matlab2tikz('width', '\fwidth', 'height', '\hwidth', '../../paper/src/figures/control-supplyair.tex',...
%     'extraaxisoptions',['xlabel style={font=\footnotesize},'...
%                        'ylabel style={font=\footnotesize},',...
%                        'legend style={font=\footnotesize},',...
%                        'ticklabel style={font=\footnotesize},'...
%                        'ylabel shift = -5 pt,'...
%                        'xlabel shift = -5 pt,']);

subplot(3,1,3); hold on; grid on; box on;
h3 = plot(t(pidx), ChwSP(pidx), '-', 'LineWidth', lw, 'Color', [0.8500    0.3250    0.0980]);

plot(52*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);
plot(56*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);
plot(64*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);
plot(68*ones(1,10), linspace(0,50,10), '--k', 'LineWidth', lw);

axis([t(pidx(1)), t(pidx(end)), 6, 9.5])
ax = gca;
ax.XTick = [40, 52, 56, 64, 68];
ax.XTickLabel = {'10am', '1pm', '2pm', '4pm', '5pm'}';
ax.YTick = [7, 9];
ylabel('temp. [^oC]')

hleg = legend(h3, 'chilled water');
set(hleg, 'box', 'on', 'Location', 'west', 'Orientation', 'Horizontal', 'EdgeColor',[1 1 1]);

% matlab2tikz('width', '\fwidth', 'height', '\hwidth', '../../paper/src/figures/control-chilled.tex',...
%     'extraaxisoptions',['xlabel style={font=\footnotesize},'...
%                        'ylabel style={font=\footnotesize},',...
%                        'legend style={font=\footnotesize},',...
%                        'ticklabel style={font=\footnotesize},'...
%                        'ylabel shift = -5 pt,'...
%                        'xlabel shift = -5 pt,']);

matlab2tikz('width', '\fwidth', 'height', '\hwidth', '../../paper/src/figures/control-all.tex',...
    'extraaxisoptions',['xlabel style={font=\footnotesize},'...
                       'ylabel style={font=\footnotesize},',...
                       'legend style={font=\footnotesize},',...
                       'ticklabel style={font=\footnotesize},'...
                       'ylabel shift = -5 pt,'...
                       'xlabel shift = -5 pt,']);
                   
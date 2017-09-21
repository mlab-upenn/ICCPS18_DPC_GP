
clear; close all;
rng(1);

%% define variables to control

SimDays = 15;
n_steps = 25;

% control variables
ClgMin = 22;
ClgMax = 32;
KitchenClgMin = 24;
KitchenClgMax = 32;
GuestClgMin = 22;
GuestClgMax = 26;
SupplyAirMin = 12;
SupplyAirMax = 14;
ChwMin = 3.7;
ChwMax = 9.7;
ctrl_vars_all = struct('ClgSP', linspace(ClgMin,ClgMax,n_steps),...
                       'KitchenClgSP', linspace(KitchenClgMin,KitchenClgMax,n_steps),...
                       'GuestClgSP', linspace(GuestClgMin,GuestClgMax,n_steps),...
                       'SupplyAirSP', linspace(SupplyAirMin,SupplyAirMax,n_steps),...
                       'ChwSP', linspace(ChwMin,ChwMax,n_steps));

% control features will be in same order
% ctrl_vars = {'ChwSP'};
ctrl_vars = {'GuestClgSP', 'SupplyAirSP', 'ChwSP'};
% ctrl_vars = {'ClgSP', 'KitchenClgSP', 'GuestClgSP', 'SupplyAirSP', 'ChwSP'};

% normalize data, except for min and max this data won't be used again
order_autoreg = 3;
n_samples = SimDays*24*4-order_autoreg;

datafile = 'unconstrained-LargeHotel';
[X, y] = load_data(datafile, order_autoreg, ctrl_vars);
[~, X_train_min, X_train_max] = preNorm(X);
[~, y_train_min, y_train_max] = preNorm(y);
offline_data = load(['../data/' datafile '.mat']);

datafile = 'test-LargeHotel';
[X_test, y_test] = load_data(datafile, order_autoreg, ctrl_vars);
X_test_norm = preNorm(X_test, X_train_min, X_train_max);
y_test_norm = preNorm(y_test, y_train_min, y_train_max);

%% create an mlepProcess instance and configure it

cd('energyPlusModels/LargeHotel/')
eplusfile = 'Building';

ep = mlepProcess;
ep.arguments = {eplusfile, 'USA_IL_Chicago-OHare.Intl.AP.725300_TMY3'};
ep.acceptTimeout = 6000; % in milliseconds

VERNUMBER = 2;  % version number of communication protocol (2 as of
                % E+ 8.1.0)

%% start EnergyPlus cosimulation
[status, msg] = ep.start;

if status ~= 0
    error('Could not start EnergyPlus: %s.', msg);
end

%% main simulation loop

EPTimeStep = 4;
deltaT = (60/EPTimeStep)*60;
kStep = 1;  % current simulation step
MAXSTEPS = SimDays*24*EPTimeStep; 

% type of random sampling
sample_type = 'prbs';
% sample_type = 'uniform';
switch sample_type
    case 'uniform'
        ClgSPrand_b7 = ClgMin+(ClgMax-ClgMin)*rand(1,MAXSTEPS);
        KitchenClgSPrand_b7 = KitchenClgMin+(KitchenClgMax-KitchenClgMin)*rand(1,MAXSTEPS);
        GuestClgSPrand_b7 = GuestClgMin+(GuestClgMax-GuestClgMin)*rand(1,MAXSTEPS);
        SupplyAirSPrand_b7 = SupplyAirMin+(SupplyAirMax-SupplyAirMin)*rand(1,MAXSTEPS);
        ChwSPrand_b7 = ChwMin+(ChwMax-ChwMin)*rand(1,MAXSTEPS);
        
        ClgSPrand_a7 = ClgMin+(ClgMax-ClgMin)*rand(1,MAXSTEPS);
        KitchenClgSPrand_a7 = KitchenClgMin+(KitchenClgMax-KitchenClgMin)*rand(1,MAXSTEPS);
        GuestClgSPrand_a7 = GuestClgMin+1+(GuestClgMax-GuestClgMin-2)*rand(1,MAXSTEPS);
        SupplyAirSPrand_a7 = SupplyAirMin+0.5+(SupplyAirMax-SupplyAirMin-1)*rand(1,MAXSTEPS);
        ChwSPrand_a7 = ChwMin+1.5+(ChwMax-ChwMin-3)*rand(1,MAXSTEPS);
        
    case 'prbs'
        ClgSPrand_b7 = postNorm(idinput(MAXSTEPS,'prbs')', ClgMin, ClgMax);
        KitchenClgSPrand_b7 = postNorm(idinput(MAXSTEPS,'prbs')', KitchenClgMin,KitchenClgMax);
        GuestClgSPrand_b7 = postNorm(idinput(MAXSTEPS,'prbs')', GuestClgMin, GuestClgMax);
        SupplyAirSPrand_b7 = postNorm(idinput(MAXSTEPS,'prbs')', SupplyAirMin, SupplyAirMax);
        ChwSPrand_b7 = postNorm(idinput(MAXSTEPS,'prbs')', ChwMin, ChwMax);
        
        ClgSPrand_a7 = postNorm(idinput(MAXSTEPS,'prbs')', ClgMin, ClgMax);
        KitchenClgSPrand_a7 = postNorm(idinput(MAXSTEPS,'prbs')', KitchenClgMin,KitchenClgMax);
        GuestClgSPrand_a7 = postNorm(idinput(MAXSTEPS,'prbs')', GuestClgMin+1, GuestClgMax-1);
        SupplyAirSPrand_a7 = postNorm(idinput(MAXSTEPS,'prbs')', SupplyAirMin+0.5, SupplyAirMax-0.5);
        ChwSPrand_a7 = postNorm(idinput(MAXSTEPS,'prbs')', ChwMin+1.5, ChwMax-1.5);
end

% variables for plotting:
outputs = nan(9,MAXSTEPS);
inputs = nan(8,MAXSTEPS);
LP = nan(1,n_samples);
RMSE = nan(1,n_samples);

% initialize: parse it to obtain building outputs
packet = ep.read;
if isempty(packet)
    error('could not read outputs from E+.');
end
[flag, eptime, outinit] = mlepDecodePacket(packet);
if flag ~= 0, error('check output flag'); end

% DOE in closed loop
iter = 0;
results = struct('chosen_x',[],'chosen_y',[]);

tic;

while kStep <= MAXSTEPS    
    if rem(kStep, 20) == 0
        fprintf('Simulation at iteration %d.\n', kStep);
    end
    
    % compute next set-points
    dayTime = mod(eptime, 86400);  % time in current day
    
    if kStep > order_autoreg
        
        iter = iter+1;
        
        % extract inputs for control features
        X_c = zeros(1,numel(ctrl_vars));
        for idc = 1:numel(ctrl_vars)
            X_c(idc) = eval(ctrl_vars{idc});
        end
        
        results.chosen_x = [results.chosen_x; ...
            preNorm([X_d, X_c], X_train_min, X_train_max)];
        results.chosen_y = [results.chosen_y; ...
            preNorm(outputs(9, kStep-1), y_train_min, y_train_max)];
        
        % training
        if rem(iter, 20) == 0 || iter==n_samples
        model = train_gp(results.chosen_x, results.chosen_y);
        results.hyperparameters(iter) = model.hyp;
        
        % save errors
        [f_star_mean_random, f_star_variance_random, ~, ~, log_probabilities] = ...
            gp(results.hyperparameters(iter), model.inference_method, ...
            model.mean_function, model.covariance_function, model.likelihood, ...
            results.chosen_x, results.chosen_y, X_test_norm, y_test_norm);
        f_star_mean_random = postNorm(f_star_mean_random, y_train_min, y_train_max);
        f_star_variance_random = postNormVar(f_star_variance_random, y_train_min, y_train_max);
        
        LP(iter) = mean(log_probabilities);
        RMSE(iter) = sqrt(mean((f_star_mean_random-y_test).^2));
        end
        
    end
    
    % disturbance features
    if kStep >= order_autoreg
        TOD = offline_data.TOD(kStep);
        DOW = offline_data.DOW(kStep);
        proxy = [TOD, DOW];
        
        Ambient = offline_data.Ambient(kStep+2-order_autoreg:kStep);
        Humidity = offline_data.Humidity(kStep-1);
        disturbances = [fliplr(Ambient'), Humidity];
        X_d = [disturbances, outputs(9, kStep-1), proxy];
        
    end
    
    % select random point
    if dayTime <= 7*3600
        ClgSP = ClgSPrand_b7(kStep);
        KitchenClgSP = KitchenClgSPrand_b7(kStep);
        GuestClgSP = GuestClgSPrand_b7(kStep);
        SupplyAirSP = SupplyAirSPrand_b7(kStep);
        ChwSP = ChwSPrand_b7(kStep);
    else
        ClgSP = ClgSPrand_a7(kStep);
        KitchenClgSP = KitchenClgSPrand_a7(kStep);
        GuestClgSP = GuestClgSPrand_a7(kStep);
        SupplyAirSP = SupplyAirSPrand_a7(kStep);
        ChwSP = ChwSPrand_a7(kStep);        
    end
    
    % let sim run for first few steps to get AR terms
    if kStep <= order_autoreg
        ClgSP = 30;
        HtgSP = 16;
        KitchenClgSP = 30;
        KitchenHtgSP = 16;
        GuestClgSP = 24;
        GuestHtgSP = 21;
        SupplyAirSP = 13;
        ChwSP = 6.7;
        
    end
    
    % need this because some inputs will follow rule-based schedules
    if dayTime <= 7*3600
        if ~any(strcmp('ClgSP',ctrl_vars)), ClgSP = 30; end;
        if ~any(strcmp('KitchenClgSP',ctrl_vars)), KitchenClgSP = 30; end;
        if ~any(strcmp('GuestClgSP',ctrl_vars)), GuestClgSP = 24; end;
        if ~any(strcmp('SupplyAirSP',ctrl_vars)), SupplyAirSP = 13; end;
        if ~any(strcmp('ChwSP',ctrl_vars)), ChwSP = 6.7; end;
        HtgSP = 16;
        KitchenHtgSP = 16;
        GuestHtgSP = 21;
        
        SP = [ClgSP, HtgSP, KitchenClgSP, KitchenHtgSP, GuestClgSP, GuestHtgSP, SupplyAirSP, ChwSP];
        
    else
        if ~any(strcmp('ClgSP',ctrl_vars)), ClgSP = 24; end;
        if ~any(strcmp('KitchenClgSP',ctrl_vars)), KitchenClgSP = 26; end;
        if ~any(strcmp('GuestClgSP',ctrl_vars)), GuestClgSP = 24; end;
        if ~any(strcmp('SupplyAirSP',ctrl_vars)), SupplyAirSP = 13; end;
        if ~any(strcmp('ChwSP',ctrl_vars)), ChwSP = 6.7; end;
        HtgSP = 21;
        KitchenHtgSP = 19;
        GuestHtgSP = 21;
        
        SP = [ClgSP, HtgSP, KitchenClgSP, KitchenHtgSP, GuestClgSP, GuestHtgSP, SupplyAirSP, ChwSP];
        
    end
    
    inputs(:,kStep) = SP;
    
    % write to inputs of E+
    ep.write(mlepEncodeRealData(VERNUMBER, 0, (kStep-1)*deltaT, SP));
    
    % read a data packet from E+
    packet = ep.read;
    if isempty(packet)
        error('Could not read outputs from E+.');
    end
    
    % initialize: parse it to obtain building outputs
    [flag, eptime, outputs(:,kStep)] = mlepDecodePacket(packet);
    if flag ~= 0, break; end
    
    kStep = kStep + 1;
    
end

% stop EnergyPlus
ep.stop;
toc;

cd('../../')

disp(['Stopped with flag ' num2str(flag)]);


%% random sampling post processing

data.Ambient = outputs(6,:)';
data.Humidity =  outputs(7,:)';
data.TotalLoad =  outputs(9,:)';
data.TOD =  outputs(3,:)';
data.DOW = outputs(4,:)';
data.ClgSP = inputs(1,:)';
data.KitchenClgSP = inputs(3,:)';
data.GuestClgSP = inputs(5,:)';
data.SupplyAirSP = inputs(7,:)';
data.ChwSP = inputs(8,:)';

[X_train, y_train] = load_data(data, order_autoreg, ctrl_vars);
X_train_norm = preNorm(X_train, X_train_min, X_train_max);
y_train_norm = preNorm(y_train, y_train_min, y_train_max);

%% results

X_chosen = X_train_norm;
y_chosen = y_train_norm;

[f_star_mean, f_star_variance, ~, ~, log_probabilities] = ...
    gp(results.hyperparameters(iter), model.inference_method, ...
       model.mean_function, model.covariance_function, model.likelihood, ...
       X_chosen, y_chosen, X_test_norm, y_test_norm);
f_star_mean = postNorm(f_star_mean, y_train_min, y_train_max);
f_star_variance = postNormVar(f_star_variance, y_train_min, y_train_max);

report = sprintf('\nRANDOM: E[log p(y* | x*, D)] = %0.3f, RMSE = %0.1f \n', ...
                 mean(log_probabilities), sqrt(mean((f_star_mean-y_test).^2)));
fprintf('%s\n', report);
loss(y_test, f_star_mean, f_star_variance);

X_chosen = postNorm(X_chosen, X_train_min, X_train_max);
y_chosen = postNorm(y_chosen, y_train_min, y_train_max);

% plotgp for random sampling
t = [0:length(y_test)-1]';
f=figure('Name', 'random sampling');
f = plotgp(f, t, y_test, f_star_mean, sqrt(f_star_variance));
axis1 = findobj(f,'Type','axes');
axis1(2).XLim = [0 size(X_test,1)];
axis1(1).XLim = [0 size(X_test,1)];

figure('Name', 'random sampling'); grid on;
yyaxis left
plot(LP, 'LineWidth', 2)
ylabel('log probability')
yyaxis right
plot(RMSE, 'LineWidth', 2)
ylabel('RMSE')

ctrl_vars_all = {'ClgSP', 'KitchenClgSP', 'GuestClgSP', 'SupplyAirSP', 'ChwSP'};
ctrl_idx = [1, 3, 5, 7, 8];
for idn = 1:numel(ctrl_vars)
    figure('Name', 'random sampling'); grid on;
    plot(inputs(ctrl_idx(strcmp(ctrl_vars{idn},ctrl_vars_all)),:), 'LineWidth', 2)
    ylabel(ctrl_vars{idn})
    xlabel('sample number')
end

%% Save results

hyperparameters = results.hyperparameters;
[YY, MM, DD, HH, MINS, ~] = datevec(now);
saveStr = sprintf('random_sampling_%s_constrained_%dinput_%dday_%04d%02d%02d_%02d%02d.mat',...
    sample_type, numel(ctrl_vars), SimDays, YY, MM, DD, HH, MINS);
save(fullfile('results', saveStr), 'model', 'hyperparameters', 'X_chosen', 'y_chosen', 'LP', 'RMSE');

save(fullfile('../data', saveStr),'-struct','data');


clear; close all;
rng(1);

[YY, MM, DD, HH, MINS, ~] = datevec(now);

%% define variables to control

SimDays = 28;       % max 28
n_steps = 20;
order_autoreg = 1+1;
n_samples = SimDays*24*4-order_autoreg;

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
ctrl_vars = {'GuestClgSP', 'SupplyAirSP', 'ChwSP'};
building = 'LargeHotel';

datafname = ['unconstrained-' building];
data_train = load(fullfile('..', 'data', datafname));
offline_data = data_train;

% Normalize the data (not all fields)
normalized_fields = {'Ambient', 'Humidity', 'TotalLoad', 'ClgSP', 'KitchenClgSP', 'GuestClgSP', 'SupplyAirSP', 'ChwSP'};
[data_train_norm, normparams] = normalize_data(data_train, normalized_fields);
y_train_min =  normparams.TotalLoad.min;
y_train_max = normparams.TotalLoad.max;

datafname = ['test-unconstrained-' building]; %'test-nominal-LargeHotel'; 'test-ramped2-LargeHotel'
data_test = load(fullfile('..', 'data', datafname));

% Normalize the data (same as for training)
data_test_norm = normalize_data(data_test, normalized_fields, normparams);

model_inputs = {...
        'TOD', ...
        'DOW', ...
        {'Ambient', order_autoreg-1:-1:0}, ...
        {'Humidity', order_autoreg-1:-1:0}, ...
        {'TotalLoad', order_autoreg-1:-1:1}, ...
        {'ClgSP', order_autoreg-1:-1:0}, ...
        {'KitchenClgSP', order_autoreg-1:-1:0}, ...
        {'GuestClgSP', order_autoreg-1:-1:0}, ...
        {'SupplyAirSP', order_autoreg-1:-1:0}, ...
        {'ChwSP', order_autoreg-1:-1:0}};
model_target = 'TotalLoad';
model_excepts = {'TOD', 'DOW'};
stepsahead = 0;

[X_train_norm, y_train_norm] = construct_data(data_train_norm, model_inputs, model_target, stepsahead, model_excepts);
[X_test_norm, y_test_norm] = construct_data(data_test_norm, model_inputs, model_target, stepsahead, model_excepts);
[~, y_test] = construct_data(data_test, model_inputs, model_target, stepsahead, model_excepts);

% search space for oed
[X_grid, Y_grid, Z_grid] = ndgrid(eval(['ctrl_vars_all.' ctrl_vars{1}]), ...
    eval(['ctrl_vars_all.' ctrl_vars{2}]), ...
    eval(['ctrl_vars_all.' ctrl_vars{3}]));
X_star = X_grid(:);
Y_star = Y_grid(:);
Z_star = Z_grid(:);
X_c_star = [X_star, Y_star, Z_star];

X_c_star_norm = zeros(size(X_c_star));
for ii = 1:numel(ctrl_vars)
    X_c_star_norm(:,ii) = preNorm(X_c_star(:,ii), normparams.(ctrl_vars{ii}).min, normparams.(ctrl_vars{ii}).max);
end
        
%% setup GP model

D = size(X_train_norm,2);

model.mean_function       = {@constant_mean};
model.covariance_function = {@ard_sqdexp_covariance};
model.likelihood          = @likGauss;

% priors on each log covariance parameter
n_cov_hyp = eval(feval(model.covariance_function{:}));
priors.cov = cell(1,n_cov_hyp);
for idc = 1:n_cov_hyp
    priors.cov{idc}  = get_prior(@gaussian_prior, 0, 1);
end

% prior on log noise
priors.lik  = {get_prior(@gaussian_prior, 0, 1)};

% prior on constant mean
priors.mean = {get_prior(@gaussian_prior, 0, 1)};

model.prior = get_prior(@independent_prior, priors);
model.inference_method = add_prior_to_inference_method(@exact_inference, model.prior);

% problem type: 'IG'-information gain or 'MV'-maximum variance
problem.type = 'IG';
problem.num_evaluations = n_samples;

%% create an mlepProcess instance and configure it

cd(['energyPlusModels/' building '/'])
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

% variables for plotting:
outputs = nan(9,MAXSTEPS);
inputs = nan(8,MAXSTEPS);
LP = zeros(1,n_samples);
RMSE = zeros(1,n_samples);
LP_map = zeros(1,n_samples);
RMSE_map = zeros(1,n_samples);

% initialize: parse it to obtain building outputs
packet = ep.read;
if isempty(packet)
    error('could not read outputs from E+.');
end
[flag, eptime, outinit] = mlepDecodePacket(packet);
if flag ~= 0, error('check output flag'); end

% DOE in closed loop
n_days = 0;
iter = 0;
tic;

while kStep <= MAXSTEPS    
    if rem(kStep, 20) == 0
        fprintf('Simulation at iteration %d.\n', kStep);
    end
    
    % compute next set-points
    dayTime = mod(eptime, 86400);  % time in current day
        
    % let sim run for first few steps to get AR terms
    if kStep <= order_autoreg
        ClgSP = 30;
        KitchenClgSP = 30;
        GuestClgSP = 24;
        SupplyAirSP = 13;
        ChwSP = 6.7;
        HtgSP = 16;
        KitchenHtgSP = 16;
        GuestHtgSP = 21;
        
        SP = [ClgSP, HtgSP, KitchenClgSP, KitchenHtgSP, GuestClgSP, GuestHtgSP, SupplyAirSP, ChwSP];
        
    else
        
        % need this because some inputs will follow rule-based schedules
        if dayTime <= 7*3600
            
            ClgSP = 30;
            KitchenClgSP = 30;
            GuestClgSP = 24;
            SupplyAirSP = 13;
            ChwSP = 6.7;
            HtgSP = 16;
            KitchenHtgSP = 16;
            GuestHtgSP = 21;
            
        else
            
            ClgSP = 24;
            KitchenClgSP = 26;
            GuestClgSP = 24;
            SupplyAirSP = 13;
            ChwSP = 6.7;
            HtgSP = 21;
            KitchenHtgSP = 19;
            GuestHtgSP = 21;
            
        end
        
        iter = iter+1;
        
        data_cur = struct();
        data_cur.TOD = offline_data.TOD(kStep);
        data_cur.DOW = offline_data.DOW(kStep);
        data_cur.Ambient = offline_data.Ambient(kStep-(order_autoreg-1):kStep)';
        data_cur.Humidity = offline_data.Humidity(kStep-(order_autoreg-1):kStep)';
        data_cur.TotalLoad = outputs(9, kStep-(order_autoreg-1):kStep-1);
        
        data_cur.ClgSP = [inputs(1,kStep-(order_autoreg-1):kStep-1), ClgSP];
        data_cur.KitchenClgSP = [inputs(3,kStep-(order_autoreg-1):kStep-1), KitchenClgSP];
        data_cur.GuestClgSP = inputs(5,kStep-(order_autoreg-1):kStep-1);
        data_cur.SupplyAirSP = inputs(7,kStep-(order_autoreg-1):kStep-1);
        data_cur.ChwSP = inputs(8,kStep-(order_autoreg-1):kStep-1);
        
        data_cur_norm = normalize_data(data_cur, normalized_fields, normparams);
        
        fields = fieldnames(data_cur_norm);
        for ff = 1:numel(fields)
            fname = fields{ff};
            data_cur_norm.(fname) = repmat(data_cur_norm.(fname), [size(X_c_star,1),1]);
            if any(strcmp(fname, ctrl_vars))
                data_cur_norm.(fname) = [data_cur_norm.(fname), X_c_star_norm(:,strcmp(fname,ctrl_vars))];
            end
        end
        
        % define candidate points for DOE (10x10x10)
        X_star_norm = [];
        for ff = 1:numel(fields)
            fname = fields{ff};
            X_star_norm = [X_star_norm, data_cur_norm.(fname)];
        end
        problem.candidate_x_star = X_star_norm;

        % select best point
        [results, ind] = learn_gp_hyperparameters_doe(problem, model, iter, results, 'num_restarts', 0);
        X_c_next = X_c_star(ind,:);
        GuestClgSP = X_c_next(1);
        SupplyAirSP  = X_c_next(2);
        ChwSP = X_c_next(3);
        
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
    
    % initial sample for experiment design
    if kStep == order_autoreg
        
        % extract features
        data_cur = struct();
        data_cur.TOD = offline_data.TOD(kStep-(order_autoreg-1):kStep);
        data_cur.DOW = offline_data.DOW(kStep-(order_autoreg-1):kStep);
        data_cur.Ambient = offline_data.Ambient(kStep-(order_autoreg-1):kStep);
        data_cur.Humidity = offline_data.Humidity(kStep-(order_autoreg-1):kStep);
        data_cur.TotalLoad = outputs(9,kStep-(order_autoreg-1):kStep);
        data_cur.ClgSP = inputs(1,kStep-(order_autoreg-1):kStep);
        data_cur.KitchenClgSP = inputs(3,kStep-(order_autoreg-1):kStep);
        data_cur.GuestClgSP = inputs(5,kStep-(order_autoreg-1):kStep);
        data_cur.SupplyAirSP = inputs(7,kStep-(order_autoreg-1):kStep);
        data_cur.ChwSP = inputs(8,kStep-(order_autoreg-1):kStep);
        
        data_cur_norm = normalize_data(data_cur, normalized_fields, normparams);
        [X_init, y_init] = construct_data(data_cur_norm, model_inputs, model_target, stepsahead, model_excepts);
        
        problem.x_init = X_init;
        problem.y_init = y_init;
        results.chosen_y = y_init;
        
    end
    
    if kStep > order_autoreg
        results.chosen_y = [results.chosen_y; ...
            preNorm(outputs(9, kStep), y_train_min, y_train_max)];
    end
    
    % save results after every 1 week
    if rem(kStep, 96*7) == 0
        n_days = n_days+7;
        
        data.Ambient = outputs(6,1:kStep)';
        data.Humidity =  outputs(7,1:kStep)';
        data.TotalLoad =  outputs(9,1:kStep)';
        data.Month = outputs(1,1:kStep)';
        data.DOM = outputs(2,1:kStep)';
        data.TOD =  outputs(3,1:kStep)';
        data.DOW = outputs(4,1:kStep)';
        data.ClgSP = inputs(1,1:kStep)';
        data.KitchenClgSP = inputs(3,1:kStep)';
        data.GuestClgSP = inputs(5,1:kStep)';
        data.SupplyAirSP = inputs(7,1:kStep)';
        data.ChwSP = inputs(8,1:kStep)';

        saveStr = sprintf('doe_%s_%s_%dinput_%dday_%04d%02d%02d_%02d%02d.mat',...
            building, problem.type, numel(ctrl_vars), n_days, YY, MM, DD, HH, MINS);
        save(fullfile('../../results', saveStr), 'model', 'results');
        save(fullfile('../../data', saveStr),'-struct','data');
        
    end
    
    kStep = kStep + 1;
    
end

% stop EnergyPlus
ep.stop;
toc;

cd('../../')

disp(['Stopped with flag ' num2str(flag)]);

%% save data

data.Ambient = outputs(6,:)';
data.Humidity =  outputs(7,:)';
data.TotalLoad =  outputs(9,:)';
data.Month = outputs(1,:)';
data.DOM = outputs(2,:)';
data.TOD =  outputs(3,:)';
data.DOW = outputs(4,:)';
data.ClgSP = inputs(1,:)';
data.KitchenClgSP = inputs(3,:)';
data.GuestClgSP = inputs(5,:)';
data.SupplyAirSP = inputs(7,:)';
data.ChwSP = inputs(8,:)';

%% DOE post processing

% errors with map estimate                                 
[f_star_mean_active, f_star_variance_active, ~, ~, log_probabilities] = ...
    gp(results.map_hyperparameters(end), model.inference_method, ...
       model.mean_function, model.covariance_function, model.likelihood, ...
       results.chosen_x, results.chosen_y, X_test_norm, y_test_norm);
   
f_star_mean_active = postNorm(f_star_mean_active, y_train_min, y_train_max);
f_star_variance_active = postNormVar(f_star_variance_active, y_train_min, y_train_max);

report_active = sprintf('\nACTIVE with MAP: E[log p(y* | x*, D)] = %0.3f, RMSE = %0.1f \n', ...
                 mean(log_probabilities), sqrt(mean((f_star_mean_active-y_test).^2)));
fprintf('%s\n', report_active);

loss(y_test, f_star_mean_active, f_star_variance_active);

% chosen samples
X_chosen = results.chosen_x;
y_chosen = results.chosen_y;

% plotgp for active learning
t = [0:length(y_test)-1]';
f=figure('Name', 'active learning');
f = plotgp(f, t, y_test, f_star_mean_active, sqrt(f_star_variance_active));
axis1 = findobj(f,'Type','axes');
axis1(2).XLim = [0 size(X_test_norm,1)];
axis1(1).XLim = [0 size(X_test_norm,1)];

figure('Name', 'active learning'); grid on;
yyaxis left
plot(LP, 'LineWidth', 2)
ylabel('log probability')
yyaxis right
plot(RMSE, 'LineWidth', 2)
ylabel('RMSE')

ctrl_vars_all = {'ClgSP', 'KitchenClgSP', 'GuestClgSP', 'SupplyAirSP', 'ChwSP'};
ctrl_idx = [1, 3, 5, 7, 8];
for idn = 1:numel(ctrl_vars)
    figure('Name', 'active learning'); grid on;
    plot(inputs(ctrl_idx(strcmp(ctrl_vars{idn},ctrl_vars_all)),:), 'LineWidth', 2)
    ylabel(ctrl_vars{idn})
    xlabel('sample number')
end

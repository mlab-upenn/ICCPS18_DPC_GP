
rng(0);

%% define variables to control

% control variables
n_steps = 15;
ctrl_variables = {'GuestClgSP', 'SupplyAirSP', 'ChwSP'};
ctrl_range = {linspace(22,26,n_steps),...
                linspace(12,14,n_steps),...
                linspace(3.7,9.7,n_steps)};

% normalize data, except for min and max this data won't be used again
datafile = 'unconstrained-LargeHotel';
order_autoreg = 3;
[X, y] = load_data(datafile, order_autoreg, ctrl_variables);
[X_norm, X_train_min, X_train_max] = preNorm(X);
[y_norm, y_train_min, y_train_max] = preNorm(y);

% offline data for future disturbances
offline_data = load(['../data/' datafile '.mat']);

[X_grid, Y_grid, Z_grid] = ndgrid(ctrl_range{1},ctrl_range{2},ctrl_range{3});
X_star = X_grid(:);
Y_star = Y_grid(:);
Z_star = Z_grid(:);
X_c_star = [X_star, Y_star, Z_star];

%% setup DOE model

model.mean_function       = {@constant_mean};
model.covariance_function = {@ard_sqdexp_covariance};
model.likelihood          = @likGauss;

% used saved initial hyperparameters
load('init_hyp.mat');

% uncomment to calculate new initial hyperparams
% n_samples_init = 1000;
% init_hyp = initial_model(file, n_samples_init, order_autoreg, ctrl_variables);

true_hyp = init_hyp;

% priors on each log covariance parameter
priors.cov = cell(1,numel(init_hyp.cov));
for idc = 1:numel(init_hyp.cov)
    priors.cov{idc}  = get_prior(@gaussian_prior, init_hyp.cov(idc), 1);
%     priors.cov{idc}  = get_prior(@gaussian_prior, 0, 1);
end

% prior on log noise
priors.lik  = {get_prior(@gaussian_prior, init_hyp.lik, 1)};
% priors.lik  = {get_prior(@gaussian_prior, 0, 1)};

% prior on constant mean
priors.mean = {get_prior(@gaussian_prior, init_hyp.mean, 1)};
% priors.mean = {get_prior(@gaussian_prior, 0, 1)};

model.prior = get_prior(@independent_prior, priors);
model.inference_method = add_prior_to_inference_method(@exact_inference, model.prior);

problem.num_evaluations = 100;

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
SimDays = 5;
deltaT = (60/EPTimeStep)*60;
kStep = 1;  % current simulation step
MAXSTEPS = SimDays*24*EPTimeStep; 

% variables for plotting:
outputs = nan(9,MAXSTEPS);
inputs = nan(8,MAXSTEPS);

% initialize: parse it to obtain building outputs
packet = ep.read;
if isempty(packet)
    error('could not read outputs from E+.');
end
[flag, eptime, outinit] = mlepDecodePacket(packet);
if flag ~= 0, error('check output flag'); end

% DOE in closed loop
iter = 0;
tic;

while kStep <= MAXSTEPS    

    % compute next set-points
    dayTime = mod(eptime, 86400);  % time in current day
    
    % disturbance features
    if kStep >= order_autoreg
        TOD = offline_data.TOD(kStep);
        DOW = offline_data.DOW(kStep);
        proxy = [TOD, DOW];
        
        Ambient = offline_data.Ambient(kStep+1-order_autoreg:kStep)';
        Humidity = offline_data.Humidity(kStep+1-order_autoreg:kStep)';
        disturbances = fliplr([Ambient; Humidity]);
        X_d = [disturbances(:); outputs(9, kStep-1); proxy']';
        
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
        
        SP = [ClgSP, HtgSP, KitchenClgSP, KitchenHtgSP, GuestClgSP, GuestHtgSP, SupplyAirSP, ChwSP];
        
    else
        
        iter = iter+1;
        
        X_d_star = repmat(X_d, [size(X_c_star,1),1]);
        
        % define candidate points for DOE (10x10x10)
        problem.candidate_x_star = [X_d_star, X_c_star];
        problem.candidate_x_star = preNorm(problem.candidate_x_star, X_train_min, X_train_max);

        % select best point
        results = learn_gp_hyperparameters_xinit(problem, model, iter, results);
        
        X_next = postNorm(results.chosen_x, X_train_min, X_train_max);
        X_c_next = X_next(end,end-2:end);
        
        GuestClgSP = X_c_next(1);
        SupplyAirSP = X_c_next(2);
        ChwSP = X_c_next(3);
        
        % need this because some inputs will follow rule-based schedules
        if dayTime <= 7*3600
            ClgSP = 30;
            HtgSP = 16;
            KitchenClgSP = 30;
            KitchenHtgSP = 16;
            % GuestClgSP = 24;
            GuestHtgSP = 21;
            % SupplyAirSP = 13;
            % ChwSP = 6.7;

            SP = [ClgSP, HtgSP, KitchenClgSP, KitchenHtgSP, GuestClgSP, GuestHtgSP, SupplyAirSP, ChwSP];
            
        else
            ClgSP = 24;
            HtgSP = 21;
            KitchenClgSP = 26;
            KitchenHtgSP = 19;
            % GuestClgSP = 24;
            GuestHtgSP = 21;
            % SupplyAirSP = 13;
            % ChwSP = 6.7;

            SP = [ClgSP, HtgSP, KitchenClgSP, KitchenHtgSP, GuestClgSP, GuestHtgSP, SupplyAirSP, ChwSP];
            
        end
        
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
        
        % extract inputs for control features
        X_c = zeros(1,numel(ctrl_variables));
        for idc = 1:numel(ctrl_variables)
            X_c(idc) = eval(ctrl_variables{idc});
        end

        X_init = [X_d, X_c];
        y_init = outputs(9, kStep);
        problem.x_init = preNorm(X_init, X_train_min, X_train_max);
        problem.y_init = preNorm(y_init, y_train_min, y_train_max);
        results.chosen_y = preNorm(y_init, y_train_min, y_train_max);
        
    end
    
    if kStep > order_autoreg
        results.chosen_y = [results.chosen_y; ...
            preNorm(outputs(9, kStep), y_train_min, y_train_max)];
    end
    
    kStep = kStep + 1;
    
end

% Stop EnergyPlus
ep.stop;
toc;

cd('../../')

disp(['Stopped with flag ' num2str(flag)]);

%% post processing

[f_star_mean_active, f_star_variance_active, ~, ~, log_probabilities] = ...
    gp(results.map_hyperparameters(end), model.inference_method, ...
       model.mean_function, model.covariance_function, model.likelihood, ...
       results.chosen_x, results.chosen_y, X_norm, y_norm);
f_star_mean_active = postNorm(f_star_mean_active, y_train_min, y_train_max);
f_star_variance_active = postNormVar(f_star_variance_active, y_train_min, y_train_max);

report_active = sprintf('ACTIVE:\n E[log p(y* | x*, D)] = %0.3f, RMSE = %0.1f', ...
                 mean(log_probabilities), sqrt(mean((f_star_mean_active-y).^2)));
fprintf('%s\n', report_active);

loss(y, f_star_mean_active, f_star_variance_active);

X_chosen_active = results.chosen_x;
y_chosen_active = results.chosen_y;
X_chosen_active = postNorm(X_chosen_active, X_train_min, X_train_max);
y_chosen_active = postNorm(y_chosen_active, y_train_min, y_train_max);

% plotgp for active learning
t = [0:length(y)-1]';
f1=figure('Name', 'active learning');
f1 = plotgp(f1, t, y, f_star_mean_active, sqrt(f_star_variance_active));
axis1 = findobj(f1,'Type','axes');
axis1(2).XLim = [0 1000];
axis1(1).XLim = [0 1000];

%% compare to random sampling


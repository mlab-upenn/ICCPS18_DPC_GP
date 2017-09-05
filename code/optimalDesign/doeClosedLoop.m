
rng(0);

%% define variables to control

% control variables
ctrl_variables = {'GuestClgSP', 'SupplyAirSP', 'ChwSP'};
ctrl_range = {linspace(22,26,10),linspace(12,14,10),linspace(3.7,9.7,10)};

% normalize data, except for min and max this data won't be used again
file = 'unconstrained-LargeHotel';
order_autoreg = 3;
[X, y] = load_data(file, order_autoreg, ctrl_variables);
[X_norm, X_train_min, X_train_max] = preNorm(X);
[y_norm, y_train_min, y_train_max] = preNorm(y);

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
% init_hyp = initial_model(file, n_samples_init);

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
        
%% create an mlepProcess instance and configure it

cd('energyPlusModels/LargeHotel/')
file = 'Building';

ep = mlepProcess;
ep.arguments = {file, 'USA_IL_Chicago-OHare.Intl.AP.725300_TMY3'};
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
SimDays = 10;
deltaT = (60/EPTimeStep)*60;
kStep = 1;  % current simulation step
MAXSTEPS = SimDays*24*EPTimeStep;  % max simulation time = 7 days

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
while kStep <= MAXSTEPS    

    % compute next set-points
    dayTime = mod(eptime, 86400);  % time in current day
    
    % disturbance features
    TOD = outputs(3, kStep-1);
    DOW = outputs(4, kStep-1);
    proxy = [TOD, DOW];
    disturbances = fliplr(outputs(1:2, kStep-order_autoreg:kStep-1));
    X_d = [disturbances(:); outputs(9, kStep-2); proxy']';
        
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
        X_d_star = repmat(Xd, [size(X_c_star,1),1]);
        problem.candidate_x_star = [X_d_star, X_c_star];
        
        results = learn_gp_hyperparameters_xinit(problem, model, iter, results);
        [GuestClgSP, SupplyAirSP, ChwSP] = results.chosen_x(end,:);
        
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
    
    % Write to inputs of E+
    ep.write(mlepEncodeRealData(VERNUMBER, 0, (kStep-1)*deltaT, SP));
    
    % Read a data packet from E+
    packet = ep.read;
    if isempty(packet)
        error('Could not read outputs from E+.');
    end
    
    % Initialize: parse it to obtain building outputs
    [flag, eptime, outputs(:,kStep)] = mlepDecodePacket(packet);
    
    % initial sample for experiment design
    if kStep == order_autoreg
        
        % extract inputs for control features
        X_c = zeros(1,numel(ctrl_variables));
        for idc = 1:numel(ctrl_variables)
            X_c(idc) = eval(ctrl_variables{idc});
        end

        X_init = [X_d, X_c];
        y_init = outputs(9, kStep);
        problem.x_init = X_init;
        problem.y_init = y_init;
        
    end
        
    if flag ~= 0, break; end

    kStep = kStep + 1;
    
end

% Stop EnergyPlus
ep.stop;

disp(['Stopped with flag ' num2str(flag)]);


figure
plot(1:MAXSTEPS,outputs(1,:));
figure
plot(1:MAXSTEPS,outputs(2,:));
figure
plot(1:MAXSTEPS,outputs(3,:));
figure
plot(1:MAXSTEPS,outputs(9,:));


% This script runs one simulation.
% Use run_experiments.m to run the actual simulations.
% Authors: Truong X. Nghiem (xuan.nghiem@epfl.ch)

%% Basic setup
% clear all
workspace = 'gpdr';
comm = 'mqtt';
mqttserver = 'tcp://localhost:1883';

% buildings = {'LargeHotel', 'MediumOffice', 'SuperMarket'};
buildings = {'LargeHotel'};
nBuildings = numel(buildings);

% The length of data used for the simulation
% datalen = (4 * 24) * (7*6);


%% Load the baseline and forecasts
BASELINE = load(fullfile('..', 'data', 'test-LargeHotel.mat'));

ref = BASELINE.TotalLoad;   % baseline in Watts
forecasts = struct('Ambient', BASELINE.Ambient, ...
    'Humidity', BASELINE.Humidity);


%% Calculate the reference tracking signal

the_baseline = ref;

% Values where hour > DRstart and <= DRend belong to the DR duration.

% Main DR event
DRidx = (BASELINE.TOD >= DRstart) & (BASELINE.TOD <= DRend);
ref(DRidx) = ref(DRidx) - DRreduction;

% after event, try to keep baseline demand
DRKickback = (BASELINE.TOD > DRend) & (BASELINE.TOD <= DRend+2);
ref(DRKickback) = ref(DRKickback);  % We may gradually ease the reference back to nominal / baseline value

% The rest = do not track
ref(~(DRidx | DRKickback)) = NaN;  % Do not track

clear BASELINE


%% Load the GP and its data

if multipleGPs
    %{
    % Load the first one
    RESULTS = load(['../' controlmodeldir '/control_cov_selection_aggregate_results_loocv.mat'], gpmodelname, 'identData');
    gpmodel = RESULTS.(gpmodelname);
    gpmodel(1).norm_ta_max = RESULTS.identData.norm_ta_max;
    gpmodel(1).norm_ta_min = RESULTS.identData.norm_ta_min;
    gpmodel(1).norm_y_min = RESULTS.identData.norm_y_min;
    gpmodel(1).norm_y_max = RESULTS.identData.norm_y_max;
    gpmodel(1).stepsahead = 0;
    %}
    
    gpmodel = struct([]);
    
    % The rest
    for k = 0:horizon-1
        gpmodelname = sprintf(['../' controlmodeldir '/gpml_final_model_ahead%02d.mat'],k);
        fprintf('Loading GP model from file: %s\n', gpmodelname);
        RESULTS = load(gpmodelname);
        gpmodel(k+1).cov = RESULTS.training_result.cov;
        if isfield(RESULTS.training_result, 'lik')
            gpmodel(k+1).lik = RESULTS.training_result.lik;
        else
            gpmodel(k+1).lik = @likGauss;
        end
        if isfield(RESULTS.training_result, 'mean')
            gpmodel(k+1).mean = RESULTS.training_result.mean;
        else
            gpmodel(k+1).mean = 'meanConst';
        end
        
        if isfield(RESULTS.training_result, 'hyp') && isstruct(RESULTS.training_result.hyp)
            gpmodel(k+1).hyp = RESULTS.training_result.hyp;
        else
            gpmodel(k+1).hyp = struct('cov', RESULTS.training_result.hypcov, ...
                'lik', log(RESULTS.training_result.gp.Sigma),...
                'mean', RESULTS.training_result.gp.Beta);
        end
        
        gpmodel(k+1).stepsahead = RESULTS.training_result.stepsahead;
        
        if isfield(RESULTS.training_result, 'gp')
            gpmodel(k+1).inputs = RESULTS.training_result.gp.X;
            gpmodel(k+1).target = RESULTS.training_result.gp.Y;
        else
            gpmodel(k+1).inputs = RESULTS.training_result.inputs;
            gpmodel(k+1).target = RESULTS.training_result.target;
        end
        
        gpmodel(k+1).lag_y = RESULTS.training_result.lag_y;
        gpmodel(k+1).lag_dr = RESULTS.training_result.lag_dr;
        gpmodel(k+1).lag_Ta = RESULTS.training_result.lag_Ta;
        gpmodel(k+1).lag_humid = RESULTS.training_result.lag_humid;

        gpmodel(k+1).norm_ta_max = RESULTS.norm_ta_max;
        gpmodel(k+1).norm_ta_min = RESULTS.norm_ta_min;
        gpmodel(k+1).norm_y_min = RESULTS.norm_y_min;
        gpmodel(k+1).norm_y_max = RESULTS.norm_y_max;
    end
    
    clear RESULTS

else
    % Load the single GP model
    RESULTS = load(fullfile('..', 'models', 'LargeHotel', control_gp_model_file));
            
    gpmodel = RESULTS.training_result.model;
    gpmodel.model_target = RESULTS.training_result.model_target;
    gpmodel.model_inputs = RESULTS.training_result.model_inputs;
    gpmodel.model_excepts = RESULTS.training_result.model_excepts;
    gpmodel.hyp = RESULTS.training_result.hyp;
    gpmodel.inputs = RESULTS.Xtrain_norm;
    gpmodel.target = RESULTS.Ytrain_norm;
    
    normparams = RESULTS.normparams;
    
    clear RESULTS
end


%% The aggregator node
% This node computes the DR signal for achieving the desired power output

% The formulation where the two-sided chance constraints are approximated
% by SOC constraints; with slack variables. A single GP is used to
% propagate to multiple steps ahead.

if ~exist('use_evolving', 'var') || ~exist('evolving_params', 'var') || ~use_evolving
    use_evolving = false;
    controller = DRtrackingWithStorage_SOCP_singlegp(gpmodel, normparams, ref, forecasts, active_ahead, horizon, ...
        wdelta, wPb, wvar,...
        ramplimit,...
        Batt_params, use_battery, ...
        'aggregator',  workspace, comm, mqttserver);
else
    use_evolving = true;
    controller = DRtrackingWithStorage_SOCP_singlegp_evolving(gpmodel, normparams, ref, forecasts, active_ahead, horizon, ...
        wdelta, wPb, wvar,...
        ramplimit,...
        Batt_params, use_battery, ...
        evolving_params, ...
        'aggregator',  workspace, comm, mqttserver);
end

disp('Please start the OBN simulation.');
try
    b = runSimulation(controller, 120, false);
catch ME
    disp('Stop simulation due to run-time error.');
    controller.stopSimulation();
    rethrow(ME);
end

%% Extract the data

if b
    pause(5);
    
    % Extract the results and save to MAT files
    RESULTS = struct('TotalLoad', 0, 'OtherLoad', 0, 'HVACLoad', 0);
    RESULTS.PMV = cell(1,numel(buildings));
    % RESULTS.CO2 = cell(1,numel(buildings));
    for kbldg = 1:numel(buildings)
        RESULTS_k = extract_results(fullfile('energyPlusModels', buildings{kbldg}, 'Output', 'Building.csv'));
        
        % Sum the power values
        RESULTS.TotalLoad = RESULTS.TotalLoad + RESULTS_k.TotalLoad;
        RESULTS.OtherLoad = RESULTS.OtherLoad + RESULTS_k.OtherLoad;
        RESULTS.HVACLoad = RESULTS.HVACLoad + RESULTS_k.HVACLoad;
        
        % Save values for weather conditions and time, only once
        if kbldg == 1
            RESULTS.Ambient = RESULTS_k.Ambient;
            RESULTS.Humidity = RESULTS_k.Humidity;
            RESULTS.Time = RESULTS_k.Time;
        end
        
        % Save PMV and CO2
        RESULTS.PMV{kbldg} = RESULTS_k.PMV;
        % RESULTS.CO2{kbldg} = RESULTS_k.CO2;
    end
    
    RESULTS.ChwSP = RESULTS_k.ChwSP;
    RESULTS.ClgSP = RESULTS_k.ClgSP;
    RESULTS.GuestClgSP = RESULTS_k.GuestClgSP;
    RESULTS.KitchenClgSP = RESULTS_k.KitchenClgSP;
    RESULTS.SupplyAirSP = RESULTS_k.SupplyAirSP;
    
    if isprop(controller, 'signals')
        signals_hist = controller.signals.getAllSignals();
        
        % The calculated setpoints (normalized)
        RESULTS.calc_ChwSP = signals_hist.ChwSP;
        RESULTS.calc_ClgSP = signals_hist.ClgSP;
        RESULTS.calc_GuestClgSP = signals_hist.GuestClgSP;
        RESULTS.calc_KitchenClgSP = signals_hist.KitchenClgSP;
        RESULTS.calc_SupplyAirSP = signals_hist.SupplyAirSP;
            
        RESULTS.PowerExpected = signals_hist.power_expected;
        RESULTS.PowerVariance = signals_hist.power_var;
        RESULTS.BatteryPower = signals_hist.batt_power;
        RESULTS.BatterySoC = signals_hist.batt_soc(2:end);
        
        if isfield(signals_hist, 'tracking_delta')
            RESULTS.TrackingDelta = signals_hist.tracking_delta;
        end
        
        if isfield(signals_hist, 'soc_confidence')
            RESULTS.SoCConfidence = signals_hist.soc_confidence;
        end
    else
        RESULTS.DRSignal = controller.dr_signal;
        RESULTS.PowerExpected = controller.power_expected;
        RESULTS.PowerVariance = controller.power_var;
        RESULTS.BatteryPower = controller.batt_power_hist;
        RESULTS.BatterySoC = controller.batt_soc_hist(2:end);
        
        if isprop(controller, 'tracking_delta')
            RESULTS.TrackingDelta = controller.tracking_delta;
        end
        
        if isprop(controller, 'soc_confidence')
            RESULTS.SoCConfidence = controller.soc_confidence;
        end
    end
    
    % Compute performance = tracking error
    result_len = length(RESULTS.TotalLoad);
    ref = ref(1:result_len);
    DRidx = DRidx(1:result_len);
    RESULTS.TrackingError = nan(result_len,1);
    RESULTS.TrackingError(DRidx) = RESULTS.TotalLoad(DRidx) + RESULTS.BatteryPower(DRidx) - ref(DRidx);
    RESULTS.ref = ref;
    RESULTS.reftol = reftol;
    
    RESULTS.active_ahead = active_ahead;
    RESULTS.horizon = horizon;
    RESULTS.wvar = wvar;
    RESULTS.wdelta = wdelta;
    RESULTS.wPb = wPb;
    RESULTS.baseline = the_baseline;
    
    RESULTS.Batt_params = Batt_params;
    RESULTS.use_battery = use_battery;
    RESULTS.DRreduction = DRreduction;
    
    if use_evolving
        RESULTS.evolving_params = evolving_params;
    end
    
    save(fullfile('results', [control_gp_model_file MATpostfix]), '-struct', 'RESULTS');
    
else
    warning('Simulation timed out!!! No data is saved.');
end

%% Plot

%{
%% Run simulations

% Run for different values of Delta and DRStart time

AllDeltas = 20:10:40;  % -20:10:40
AllStartTimes = 8:14;

% Requested change in power demand, in kW
%Delta = 20;
%DRstart = 12;    % the hour at which DR starts (hour), everyday

for Delta = AllDeltas
    for DRstart = AllStartTimes
        fprintf('Simulation for Delta = %g, start time = %g\n', Delta, DRstart);
        
        DRduration = 4; % duration of DR, in hours
        DRend = DRstart + DRduration;

        aggregatornode = DRtracking(GPDATA.covFunc, GPDATA.hyp, GPDATA.lik, autoregressive_lengths,...
            GPDATA.xtrain, GPDATA.ytrain, GPDATA.avgy,...
            DRstart, DRend, baseline, Delta,...
            'aggregator', workspace, 'mqtt', mqttserver);
        

        % The building node
        drctrlnodes = createDRNode('buildingctrl', workspace, 'mqtt', mqttserver);


        % Run simulation
        %fprintf('Starting the simulation: please start all other nodes then the SMN...');
        
        % For debugging purposes, we shall run the simulation with a large timeout
        % value (so that if some other node hangs, this one will eventually time
        % out and will not hang. We will also turn off stopping the simulation if
        % timeout, so that later on we can resume the simulation if desired.
        %profile on -timer cpu
        b = runSimulation([aggregatornode drctrlnodes], 30, false);
        %profile viewer
        %fprintf('Simulation stopped with result (1/true means success): %d.\n', b);
        
        
        % Extract results and save
        result_data = struct('drsignal', aggregatornode.dr_signal,...
            'power_expected', aggregatornode.power_expected,...
            'power_var', aggregatornode.power_var,...
            'baseline', baseline,...
            'delta', Delta, 'starttime', DRstart);
        
        % pause to make sure the E+ file was written
        pause(2);
        ktrial = 0;
        while ktrial < 10
            if exist('Output/MediumOffice.csv', 'file') > 0
                break;
            else
                pause(1);
            end
            ktrial = ktrial + 1;
        end
        
        [~,result_data.power_actual,~] = mlepLoadEPResults('Output/MediumOffice.csv', 8);
        result_data.power_actual = result_data.power_actual * 1e-3;
        result_data.delta_actual = result_data.power_actual - baseline(1:length(result_data.power_actual));
        
        result_file_name = sprintf('result_delta%g_start%g', Delta, DRstart);
        save(result_file_name, '-struct', 'result_data');
        
        % Clear the nodes to start again
        clear aggregatornode drctrlnodes
        pause(2);
    end
end

% figure; hold on;
% n = length(aggregatornode.power_expected) - 1;
% hl1 = boundedline(1:n, aggregatornode.power_expected(1:n) + GPDATA.avgy, 2*sqrt(aggregatornode.power_var(1:n)), 'k');
% hl2 = plot(baseline(1:n));
% [VARS,DATA,TS] = mlepLoadEPResults('Output/MediumOffice.csv', [8]);
% hl3 = plot(DATA*1e-3, 'r');
% xlabel('Time step');
% ylabel('Power [kW]');
% legend([hl1 hl2 hl3], 'Expected', 'Baseline', 'Actual');
% 
% figure; hold on;
% plot(DATA*1e-3 - baseline(1:n));
% plot([1 n], [Delta Delta], '--');
% legend('Actual', 'Requested');
% title('Actual vs requested power change');

% delete the node variable to close it
% clear node
%}
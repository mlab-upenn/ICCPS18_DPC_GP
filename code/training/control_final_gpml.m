% Train final GPs with selected features using GPML v4.0
%

%% Load and prepare data
RESULTS = load('../Aggregate_Chicago_data_chiller_iidunif');
matprefix = 'gpml_';
matinfolder = 'aggregate_chiller_iidunif_2week';

dr_timestep = 15;       % The time step in minutes
remove_nonworkdays = true;     % If we want to remove non-work days
                               % (Sundays and holidays)
%working_hours = [8,22];     % working hours, during which we want to learn the models
working_hours = [];

stepsperhour = 60 / dr_timestep;

% Construct the identification data and validation data (for feature
% selection).

% Number of data points to use for identification
ident_len = 2*7 * (24*stepsperhour);

% Number of data points to use for validation (right after the ident data)
val_len = 1*7 * (24*stepsperhour);

% Some of the data values are shifted by one step (by EnergyPlus), e.g.
% Power[k+1] actually corresponds to the time step t_k (hence u[k]).
% However, as GPDyn (specifically its construct() function) assumes the
% model: y[k] = f(y[k-1], ... , y[k-lag_y], u[k-1], ... , u[k-lag_u])
% we should not shift these data values when loading them from the MAT file

% Obtain the ident data
identData = struct();

identData.u_hour = RESULTS.Time(1:ident_len, 3:4) * [1; 1/60];  % 0-24
% The day of week is shifted by one and is an input
identData.u_day = RESULTS.Time((1:ident_len)+1, 6);  % 1-7
if remove_nonworkdays
    % not used in GP; also shifted by one-step
    identData.u_holiday = RESULTS.Time((1:ident_len)+1, 7);
end

identData.u_ta = RESULTS.Ambient(1:ident_len,:);
identData.u_humid = RESULTS.Humidity(1:ident_len,:) / 100;  % 0 - 1

% DR signal is shifted 1 step
identData.u_dr = RESULTS.DRSignal((1:ident_len)+1,:);  % -1 - 1

% Output: interval demand, it is shifted one time step by EnergyPlus, but
% we must load it as-is (see note above).
identData.y = RESULTS.TotalLoad(1:ident_len,:) * 1e-3;  % in kW

% Obtain the validation data
% See notes above regarding indices, etc.
valData = struct();

valData.u_hour = RESULTS.Time((1:val_len)+ident_len, 3:4) * [1; 1/60];  % 0-24
% The day of week is shifted by one and is an input
valData.u_day = RESULTS.Time((1:val_len)+1+ident_len, 6);  % 1-7
if remove_nonworkdays
    % not used in GP
    valData.u_holiday = RESULTS.Time((1:val_len)+1+ident_len, 7);
end

valData.u_ta = RESULTS.Ambient((1:val_len)+ident_len,:);
valData.u_humid = RESULTS.Humidity((1:val_len)+ident_len,:) / 100;  % 0 - 1

% DR signal is shifted 1 step
valData.u_dr = RESULTS.DRSignal((1:val_len)+1+ident_len,:);  % -1 - 1

valData.y = RESULTS.TotalLoad((1:val_len)+ident_len,:) * 1e-3;  % in kW


% Normalizing data
% We don't normalize time inputs (hour and day) because they are small and
% will enter periodic covariance functions (which implicitly convert the
% input to [-1,1] range.
%
% We only normalize inputs that have large ranges.
%
% For the output values, we need to normalize them before constructing the
% training inputs (which include lagged output values) and targets, because
% the normalization mapping of the output must be the same for consistent
% multi-step simulation (it doesn't matter for 1-step prediction).

% Normalize identification data and calculate the min and max values
[identData.u_ta_norm, norm_ta_min, norm_ta_max] = preNorm(identData.u_ta);
[identData.y_norm, norm_y_min, norm_y_max] = preNorm(identData.y);

identData.norm_ta_min = norm_ta_min;
identData.norm_ta_max = norm_ta_max;
identData.norm_y_min = norm_y_min;
identData.norm_y_max = norm_y_max;

% Normalize the validation data using the computed min and max values
valData.u_ta_norm = preNorm(valData.u_ta, norm_ta_min, norm_ta_max);
valData.y_norm = preNorm(valData.y, norm_y_min, norm_y_max);  % Used to construct AR inputs for 1-step predictions


% Inline if
iif = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();

remove_humidity = false;

%% Train and validate for each stepsahead value

for stepsahead = 0:7
    matinfilename = sprintf('%s/gpml_feature_selection_ahead%02d', matinfolder, stepsahead);
    matoutfilename = sprintf('%s/gpml_final_model_ahead%02d', matinfolder, stepsahead);
    
    % Load the feature selection file to get configuration
    load(matinfilename);
    lag_y = selection_result.lag_y;
    lag_dr = selection_result.lag_dr;
    lag_Ta = selection_result.lag_Ta;
    lag_humid = selection_result.lag_humid;
        
    Nins = numel(lag_y) + numel(lag_Ta) + numel(lag_humid) + numel(lag_dr);
    
    %{
    % SE + RQ
    % Construct the covariance function
    cov = {'covProd', {...
        {'covMask', {1:Nins, 'covSEard'}},...
        {'covMask', {Nins+1, 'covRQiso'}}} ...  % temporal
        };
    
    hypnames = [arrayfun(@(d) ['y', num2str(d)], lag_y, 'UniformOutput', false),...
        arrayfun(@(d) ['u', num2str(d)], lag_dr, 'UniformOutput', false),...
        arrayfun(@(d) ['Ta', num2str(d)], lag_Ta, 'UniformOutput', false),...
        arrayfun(@(d) ['h', num2str(d)], lag_humid, 'UniformOutput', false),...
        {'sf', 't', 'sf', 'alpha'}]';
    %}

    % (SE + Const)*RQ
    % Construct the covariance function
    cov = {'covProd', {...
        {'covMask', {1:Nins, {'covSum', {'covConst', 'covSEard'}}}},...
        {'covMask', {Nins+1, 'covRQiso'}}} ...  % temporal
        };
    
    hypnames = [{'const'},...
        arrayfun(@(d) ['y', num2str(d)], lag_y, 'UniformOutput', false),...
        arrayfun(@(d) ['u', num2str(d)], lag_dr, 'UniformOutput', false),...
        arrayfun(@(d) ['Ta', num2str(d)], lag_Ta, 'UniformOutput', false),...
        arrayfun(@(d) ['h', num2str(d)], lag_humid, 'UniformOutput', false),...
        {'sf', 't', 'sf', 'alpha'}]';
    
%{    
    % (SE + Const)*RQ(periodic)
    % Construct the covariance function
    cov = {'covProd', {...
        {'covMask', {1:Nins, {'covSum', {'covConst', 'covSEard'}}}},...
        {'covMask', {Nins+1, {'covPERiso', {'covRQiso'}}}}} ...  % temporal
        };
    
    hypnames = [{'const'},...
        arrayfun(@(d) ['y', num2str(d)], lag_y, 'UniformOutput', false),...
        arrayfun(@(d) ['u', num2str(d)], lag_dr, 'UniformOutput', false),...
        arrayfun(@(d) ['Ta', num2str(d)], lag_Ta, 'UniformOutput', false),...
        arrayfun(@(d) ['h', num2str(d)], lag_humid, 'UniformOutput', false),...
        {'sf', 'per', 't', 'sf', 'alpha'}]';
        
    % For periodic kernel for temporal input (time of day), we need an
    % additional hyperparameter 'per', which does not exist in the original
    % feature selection result (as of this coding).
    if isstruct(selection_result.hyp)
        my_hyps = selection_result.hyp;
        my_hyps.cov = [vec(selection_result.hyp.cov(1:end-3)); log(24); vec(selection_result.hyp.cov(end-2:end))];
    elseif isvector(selection_result.hyp)
        my_hyps = [vec(selection_result.hyp(1:end-3)); log(24); vec(selection_result.hyp(end-2:end))];
    else
        error('Unknown hyperparameter data type.');
    end
%}
        
    my_hyps = selection_result.hyp;
        
    % Note that we initialize the hyperparameters to the ones in feature
    % selection
    training_result = control_train_gpml(identData, lag_y, lag_Ta, lag_humid, lag_dr, cov, remove_nonworkdays, my_hyps, stepsahead, working_hours);
    
    training_result.lag_y = lag_y;
    training_result.lag_Ta = lag_Ta;
    training_result.lag_humid = lag_humid;
    training_result.lag_dr = lag_dr;
    training_result.stepsahead = stepsahead;
    training_result.hypnames = hypnames;
    training_result.hours = working_hours;
    
    % Validate by performing 1-step predictions, then compute several measures
    % compared to the real outputs
    validation_result = control_validation_gpml(training_result,valData,lag_y,lag_Ta,lag_humid,lag_dr,remove_nonworkdays,stepsahead,norm_y_min,norm_y_max, working_hours);
    
    save(matoutfilename, 'training_result', 'validation_result', 'norm_ta_min', 'norm_ta_max', 'norm_y_min', 'norm_y_max');
    
end
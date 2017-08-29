function [inputs, target] = construct_data(data, lag_y, lag_Ta, lag_humid, lag_dr, remove_nonworkdays, stepsahead, testtarget, hours)
%CONSTRUCT_DATA Construct input, output data for control GP model.
%
%Inputs:
%   data : structure containing the source data fields
%   testtarget: optional vector of target outputs, in case it's not in data
%
% The inputs are ordered as followed:
% [ar_outputs, dr_signals, Ta_inputs, humid_inputs, temporal_inputs]
% where
%   dr_signals = [s(k-L1), s(k-L12),..., s(k-L1p)] where Li in lag_dr
%   ar_outpus = [y(k-L1), y(k-L2),..., y(k-Lp)] where Li in lag_y
%   Ta_inputs = [Ta(k-L1), Ta(k-L2),...,Ta(k-Lp)] where Li in lag_Ta
%   humid_inputs = [h(k-L1), h(k-L2),...,h(k-Lp)] where Li in lag_humid
%   temporal_inputs = [hour_of_day, day_of_week]
%
% If stepsahead is provided (default = 0), then the output of the model is
% for stepsahead time steps ahead of the current time, all lag values are
% counted from (t + stepsahead), however the time-of-day input to the GP is
% t, not (t + stepsahead). For example if current time is t = 2.0 and with
% sampling time of 15 minutes, stepsahead = 2, then all lag values are
% counted from 2h30' while the time input to the GP is t=2, not t = 2.5.
%
% hours, if present, will specify the hours we take into the learning; it's
% of the form [starthour, endhour] for starthour <= t < endhour.


has_testtarget = exist('testtarget', 'var') && isvector(testtarget);
if has_testtarget
    assert(numel(testtarget) == numel(data.y_norm),...
        'Test target must have the same length as of data.');
end

if exist('stepsahead', 'var') && isscalar(stepsahead) && isnumeric(stepsahead)
    assert(stepsahead >= 0, 'Number of steps ahead must be non-negative.');
else
    stepsahead = 0;
end

% Process lag values
lag_y = unique(lag_y);
assert(~isempty(lag_y));
assert(min(lag_y) >= 1);
max_lag_y = max(lag_y);

lag_Ta = unique(lag_Ta);
assert(~isempty(lag_Ta));
assert(min(lag_Ta) >= 1);
max_lag_Ta = max(lag_Ta);
excluded_Ta = ~ismember(max_lag_Ta:-1:1, lag_Ta);

if ~isempty(lag_humid)
    lag_humid = unique(lag_humid);
    assert(~isempty(lag_humid));
    assert(min(lag_humid) >= 1);
    max_lag_humid = max(lag_humid);
    excluded_humid = ~ismember(max_lag_humid:-1:1, lag_humid);
else
    max_lag_humid = 0;
    excluded_humid = true(0);
end

lag_dr = unique(lag_dr);
assert(~isempty(lag_dr));
assert(min(lag_dr) >= 1);
max_lag_dr = max(lag_dr);
excluded_dr = ~ismember(max_lag_dr:-1:1, lag_dr);

max_lag = max([max_lag_y, max_lag_Ta, max_lag_humid, max_lag_dr, stepsahead]);
excluded_y = ~ismember(max_lag:-1:1, lag_y);

% Indices that will be excluded
excluded_inputs = [excluded_y(:); excluded_dr(:); excluded_Ta(:); excluded_humid(:)];
included_inputs = find(~excluded_inputs);

% Construct the inputs based on the lag values, for DR inputs
[inputs, target] = construct([max_lag, max_lag_dr], data.u_dr, data.y_norm);

% Do it separately for the Ta input; append it to the other inputs.
Ta_inputs = construct([max_lag, max_lag_Ta], data.u_ta_norm);
inputs = [inputs, Ta_inputs];

% Do it separately for the humidity input; append it to the other
% inputs.
if ~isempty(lag_humid)
    humidinputs = construct([max_lag, max_lag_humid], data.u_humid);  % Flip the lag values to get the correct output matrix for humid inputs
    inputs = [inputs, humidinputs];
end

% The final order of inputs is: AR-power-outputs, DR-inputs, 
% Ta-inputs, (optional) humid-inputs, temporal-inputs

% Remove those that we don't need
inputs = inputs(:, included_inputs);    % Not the most efficient way, but works for CasADi

% nLaggedInputs = size(inputs, 2);
nPoints = size(inputs, 1);   % number of training data points

% Concatenate with the temporal inputs at the current time t, not
% (t+stepsahead)
inputs = [inputs, data.u_hour(end-nPoints-stepsahead:end-1-stepsahead,:)]; %, data.u_day(end-nPoints:end-1,:)];

if has_testtarget
    target = testtarget(end-nPoints+1:end);
end

% Remove non-work days' data if desired
if remove_nonworkdays
    kept_idx = find(~((data.u_day(end-nPoints:end-1,end) < 2) | (data.u_day(end-nPoints:end-1,end) > 6) | (data.u_holiday(end-nPoints:end-1,:) ~= 0)));
    inputs = inputs(kept_idx,:);    % Not the most efficient way, but works for CasADi
    target = target(kept_idx,:);
end

if exist('hours', 'var') && isvector(hours) && length(hours) == 2
    starthour = hours(1);
    endhour = hours(2);
    assert(starthour < endhour, 'Invalid hours provided.');
    
    % Remove those outside [starthour, endhour)
    kept_idx = inputs(:,end) >= starthour & inputs(:,end) < endhour;
    inputs = inputs(kept_idx,:);
    target = target(kept_idx,:);
end

end

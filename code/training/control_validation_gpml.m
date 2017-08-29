function valResults = control_validation_gpml(trainResults, valData, lag_y, lag_Ta, lag_humid, lag_dr, remove_nonworkdays, stepsahead, norm_y_min, norm_y_max, hours)
%trainResults is the result structure returned by the training function
%control_*train.m
% lag_* are the lagged outputs and non-temporal inputs.
% Each is a vector of lag values, and only those are considered as inputs.
% 
% The inputs are ordered as followed:
% [ar_outputs, dr_signals, Ta_inputs, humid_inputs, temporal_inputs]
% where
%   dr_signals = [s(k-L1), s(k-L12),..., s(k-L1p)] where Li in lag_dr
%   ar_outpus = [y(k-L1), y(k-L2),..., y(k-Lp)] where Li in lag_y
%   Ta_inputs = [Ta(k-L1), Ta(k-L2),...,Ta(k-Lp)] where Li in lag_Ta
%   Ta_inputs = [h(k-L1), h(k-L2),...,h(k-Lp)] where Li in lag_humid
%   temporal_inputs = [hour_of_day, day_of_week]
%
% For stepsahead, see construct_data() in this folder.
%
% hours, if present, will specify the hours we take into the learning; it's
% of the form [starthour, endhour] for starthour <= t < endhour.


if ~exist('stepsahead', 'var')
    stepsahead = 0;
end

valResults = struct();
if ~exist('hours', 'var')
    hours = [];
end
[valResults.inputs, valResults.target] = construct_data(valData, lag_y, lag_Ta, lag_humid, lag_dr, remove_nonworkdays, stepsahead, valData.y, hours);

[valResults.y_norm, valResults.s2_norm] = gp(trainResults.hyp, @infExact, trainResults.mean, trainResults.cov, trainResults.lik,...
    trainResults.inputs, trainResults.target, valResults.inputs);

% De-normalize the data
valResults.y = postNorm(valResults.y_norm, norm_y_min, norm_y_max);
valResults.s2 = postNormVar(valResults.s2_norm, norm_y_min, norm_y_max);

% Compute some loss measures
[valResults.loss_ae, valResults.loss_se, valResults.loss_lpd, valResults.loss_mrse, valResults.loss_smse, valResults.loss_msll] = ...
    loss(valResults.target, valResults.y, valResults.s2);

end

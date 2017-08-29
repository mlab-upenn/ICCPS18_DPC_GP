function results = control_train_gpml(identData, lag_y, lag_Ta, lag_humid, lag_dr,...
    cov, remove_nonworkdays, mlhyps, stepsahead, hours)
    
% Train an ARX GP model for control.
% covfunc is the full covariance function for all inputs, using GPML.
% The inputs are ordered as followed:
% [ar_outputs, dr_signals, Ta_inputs, humid_inputs, temporal_inputs]
% where
%   dr_signals = [s(k-L1), s(k-L12),..., s(k-L1p)] where Li in lag_dr
%   ar_outpus = [y(k-L1), y(k-L2),..., y(k-Lp)] where Li in lag_y
%   Ta_inputs = [Ta(k-L1), Ta(k-L2),...,Ta(k-Lp)] where Li in lag_Ta
%   Ta_inputs = [h(k-L1), h(k-L2),...,h(k-Lp)] where Li in lag_humid
%   temporal_inputs = [hour_of_day, day_of_week]
%
% For stepsahead, see construct_data.m in this folder.
%
% hours, if present, will specify the hours we take into the learning; it's
% of the form [starthour, endhour] for starthour <= t < endhour.

if ~exist('stepsahead', 'var')
    stepsahead = 0;
end

% Construct the GP model's components
if isstruct(mlhyps)
    % We can provide the full hyp structure
    hyp = mlhyps;
elseif isvector(mlhyps)
    hyp = struct();
    hyp.cov = mlhyps;
else
    error('Invalid hyp value: must be either a vector or a struct.');
end

% Likelihood (output noise)
lik = 'likGauss';
if ~isfield(hyp, 'lik')
    hyp.lik = log(0.1);
end

% Mean function
meanfunc = 'meanConst';
if ~isfield(hyp, 'mean')
    hyp.mean = 0;
end

theinf = @infExact;

% Construct the data
if ~exist('hours', 'var')
    hours = [];
end
[inputs, target] = construct_data(identData, lag_y, lag_Ta, lag_humid, lag_dr, remove_nonworkdays, stepsahead, [], hours);

% Validate the provided hyperparamers
D = size(inputs, 2);
validate_gp(hyp, theinf, meanfunc, cov, lik, D);

% Train the model
[hyp, flogtheta, ~] = trainGParx(hyp, theinf, meanfunc, cov, lik, inputs, target, @minimize_minfunc, -100);

% Return the results
results = struct();
results.hyp = hyp;
results.cov = cov;
results.lik = lik;
results.mean = meanfunc;
results.inputs = inputs;
results.target = target;
results.flogtheta = flogtheta;

end

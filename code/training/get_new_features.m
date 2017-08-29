function [lag_y, lag_dr, lag_Ta, lag_humid, hyp] = get_new_features(training_result, cutoff)
% Returns the features that have log lengthscales <= cutoff value.
% Also returns new lags.

nhyps = length(training_result.hyp.cov);

len_y = length(training_result.lag_y);
len_dr = length(training_result.lag_dr);
len_Ta = length(training_result.lag_Ta);
len_humid = length(training_result.lag_humid);

Nins = len_y + len_dr + len_Ta + len_humid;

if nhyps > Nins+5
    % temporal kernel is periodic
    nhyplast = 5;
else
    % temporal kernel is not periodic
    nhyplast = 4;
end
idx = (2:nhyps-nhyplast); % indices of lengthscales

ls = training_result.hyp.cov(idx);

keep = ls <= cutoff;    % which hyps we want to keep

idxkept = idx(keep);    % indices in the original hyp vector we want to keep

% Print them out
disp('Remaining features:');
fprintf('%s,', training_result.hypnames{idxkept});
fprintf('\n');

% Calculate new lags
lag_y = training_result.lag_y(keep(1:len_y));
keep(1:len_y) = []; % remove them

lag_dr = training_result.lag_dr(keep(1:len_dr));
keep(1:len_dr) = []; % remove them

lag_Ta = training_result.lag_Ta(keep(1:len_Ta));
keep(1:len_Ta) = []; % remove them

lag_humid = training_result.lag_humid(keep(1:len_humid));
keep(1:len_humid) = []; % remove them

assert(isempty(keep));

% Get new hyp values (can be used to initialize the new GP)
hyp = training_result.hyp;
hyp.cov = hyp.cov([1, idxkept, (end-nhyplast+1):end]);

end
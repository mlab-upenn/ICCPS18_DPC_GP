function select_features(matfile)
% An interactive function to select features from a result file (created
% with GPML).

assert(ischar(matfile));

load(matfile, 'training_result');

Nins = numel(training_result.lag_y) + numel(training_result.lag_Ta) + numel(training_result.lag_humid) + numel(training_result.lag_dr);

% Plot lengthscales, any values larger than 10 will be set to 10 to avoid
% too large values
disp('Plotting the log lengthscales');
f = figure;

if length(training_result.hyp.cov) > Nins + 5
    % temporal kernel is periodic
    hyps = training_result.hyp.cov(2:end-5);
else
    % temporal kernel is not periodic
    hyps = training_result.hyp.cov(2:end-4);
end
hyps(hyps > 10) = 10;
barh(hyps);

cutoff = input('Please enter the cutoff value: ');
close(f);
if isempty(cutoff) || ~isnumeric(cutoff)
    error('A valid cutoff value must be entered.');
end

% Select the features and save
selection_result = struct();
[selection_result.lag_y, selection_result.lag_dr, selection_result.lag_Ta, selection_result.lag_humid, selection_result.hyp] = get_new_features(training_result, cutoff);
save(matfile, 'selection_result', '-append');

end
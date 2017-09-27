function [X, Y] = construct_data(data, inputs, target, stepsahead, excepts)
%CONSTRUCT_DATA Construct input, output data for control GP model.
%
%Inputs:
%   data : structure containing the source data fields, each field is a
%           column vector
%   inputs : a cell array that specifies the inputs / features. Each cell
%           is either:
%           + a cell {name, lags} where name is the field of the input in data
%               structure (must be a valid field name), and lags is a
%               vector of lag values as in lagmatrix.m. See the notes below.
%           + a single string which is equivalent to {name, 0}
%           
%           The order of the inputs will be kept as is.
%   target : the field name of the target / output, in the data structure.
%
% If stepsahead is provided (default = 0), then the output of the model is
% for stepsahead time steps ahead of the current time, all lag values are
% counted from (t + stepsahead), however see the exception below. For
% example if current time is t = 2.0 and with sampling time of 15 minutes,
% stepsahead = 2, then all lag values are counted from 2h30', i.e., lag
% value of 0 is at 2h30', lag value of 1 is at 2h15', and so on.
%
% excepts (default = {}) is a cell array of names of inputs that are not
% subject to the time shift due to stepsahead mentioned above. Using the
% same example above but excepts = {'TOD'} (time of day) then for input
% 'TOD', lag value of 0 is at 2h, not 2h30'.
%
% X and Y are the constructed matrices of inputs and output for the model.
% If there is any non-zero lag value, those rows with NaN (missing values)
% will be deleted from X and Y, so X and Y might be shorter than the
% original data.

%% Check arguments
assert(isstruct(data));
allinputs = fieldnames(data);
assert(~isempty(allinputs), 'DATA must contain at least one field.');

assert(all(structfun(@isvector, data)), 'All fields of DATA must be vectors.');
datalen = unique(structfun(@length, data));
assert(isscalar(datalen), 'All fields of DATA must have the same length.');
datalen = datalen(1);

assert(iscell(inputs), 'INPUTS must be a cell array.');
% Make inputs have a standard form {{name, lags}, ...}; also check the
% names
totalinputs = 0;
for k = 1:numel(inputs)
    if ischar(inputs{k})
        assert(isfield(data, inputs{k}), 'The input "%s" is not a valid field in DATA.', inputs{k});
        inputs{k} = {inputs{k}, 0};
        totalinputs = totalinputs + 1;
    else
        assert(iscell(inputs{k}) && numel(inputs{k}) == 2 && ...
            ischar(inputs{k}{1}) && isnumeric(inputs{k}{2}) && ...
            isvector(inputs{k}{2}), ...
            'Invalid element %d of INPUTS.', k);
        assert(isfield(data, inputs{k}{1}), 'The input %d with name "%s" is not a valid field in DATA.', k, inputs{k}{1});
        assert(~isempty(inputs{k}{2}) && all(inputs{k}{2} >= 0), 'Invalid lags for input %d.', k);
        totalinputs = totalinputs + numel(inputs{k}{2});
    end
end

if exist('stepsahead', 'var') && isscalar(stepsahead) && isnumeric(stepsahead)
    assert(stepsahead >= 0, 'Number of steps ahead must be non-negative.');
else
    stepsahead = 0;
end

if exist('excepts', 'var')
    if ischar(excepts)
        excepts = {excepts};
    else
        assert(iscellstr(excepts), 'EXCEPTS must be a string or a cell array of string.');
    end
else
    excepts = {};
end
if ~isempty(excepts)
    assert(all(ismember(excepts, allinputs)), 'EXCEPTS must only contain fields in DATA.');
end

assert(ischar(target) && ismember(target, allinputs), 'TARGET must be the name of a field in DATA.');

%% Loop through the inputs and construct the input matrix X

X = zeros(datalen, totalinputs);    % pre-allocate X matrix
curIdx = 1;
for k = 1:numel(inputs)
    thelags = inputs{k}{2};
    numlags = numel(thelags);
    thefield = inputs{k}{1};
    if ~ismember(thefield, excepts)
        thelags = thelags - stepsahead;
    end
    X(:, curIdx:curIdx+numlags-1) = lagmatrix(vec(data.(thefield)), thelags);
    curIdx = curIdx + numlags;
end
assert(curIdx == totalinputs + 1);

%% Construct the output Y
Y = lagmatrix(vec(data.(target)), -stepsahead);

%% Remove NaN samples
Xnan = any(isnan(X), 2);    % Xnan = column vector indicating NaN rows in X
nanRows = Xnan | isnan(Y);  % NaN rows in X or Y

% Remove these rows
X(nanRows, :) = [];
Y(nanRows) = [];

end

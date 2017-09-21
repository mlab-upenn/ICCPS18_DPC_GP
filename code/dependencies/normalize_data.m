function [data, normparams] = normalize_data(data, fields, normparams)
%NORMALIZE_DATA Normalize given fields of a data structure.
%
% INPUTS:
%   data : data structure.
%   fields : optional list of fields to be normalized.
%   normparams : optional normalization parameter structure, as returned by
%       this function. If it is provided, the data will be normalized
%       according to these parameters.
%
% OUTPUTS:
%   datanorm : new data structure with given fields normalized.
%   normparams : structure of the min and max values used for normalizing
%       the given fields.

assert(isstruct(data));

if nargin < 2 || isempty(fields), fields = fieldnames(data); end
if ischar(fields), fields = {fields}; end
assert(iscellstr(fields));

if nargin < 3 || isempty(normparams), normparams = struct; end
assert(isstruct(normparams));

for k = 1:numel(fields)
    fname = fields{k};
    if isfield(data, fname)
        % If the normalization parameters are given, use them
        if isfield(normparams, fname)
            data.(fname) = preNorm(data.(fname), normparams.(fname).min, normparams.(fname).max);
        else
            % Else, calculate them
            [data.(fname), normparams.(fname).min, normparams.(fname).max] = preNorm(data.(fname));
        end
    end
end

end


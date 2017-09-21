function data = denormalize_data(data, normparams, fields)
%DENORMALIZE_DATA De-normalize a data structure.
%
% INPUTS:
%   data : data structure.
%   normparams : the normalization parameter structure, as returned by
%               NORMALIZE_DATA.
%   fields : optional list of fields to be de-normalized.
%
% OUTPUTS:
%   data : new data structure with de-normalized fields.

assert(isstruct(data));
assert(isstruct(normparams));

if nargin < 3 || isempty(fields), fields = fieldnames(normparams); end
if ischar(fields), fields = {fields}; end
assert(iscellstr(fields));

for k = 1:numel(fields)
    fname = fields{k};
    if isfield(normparams, fname) && isfield(data, fname)
        data.(fname) = postNorm(data.(fname), normparams.(fname).min, normparams.(fname).max);
    end
end

end
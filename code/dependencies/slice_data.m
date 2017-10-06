function data = slice_data(data, len)
% Slide fields of a data structure.
fldnames = fieldnames(data);
toend = isscalar(len);
for k = 1:numel(fldnames)
    if toend
        data.(fldnames{k}) = data.(fldnames{k})(len:end,:);
    else
        data.(fldnames{k}) = data.(fldnames{k})(len(1):len(2),:);
    end
end
end
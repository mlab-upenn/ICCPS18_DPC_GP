function data = slice_data(data, len)
% Slide fields of a data structure.
fldnames = fieldnames(data);
for k = 1:numel(fldnames)
    data.(fldnames{k}) = data.(fldnames{k})(len(1):len(2),:);
end
end
function data = slice_data_batch(data, idx)
% Slide fields of a data structure.
fldnames = fieldnames(data);
for k = 1:numel(fldnames)
    data.(fldnames{k}) = data.(fldnames{k})(idx,:);
end
end
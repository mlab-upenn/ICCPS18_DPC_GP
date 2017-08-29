function HV = get_hyp_values(RESULT, HYP)
% Returns the values of hyperparameters in RESULT whose names start with
% HYP.
% HYP is a string, empty to get all hyperparameters.
% RESULT is a structure with field 'hypcov' and 'hypnames'.

HV = struct('name', {}, 'value', {});
getAll = isempty(HYP);
if getAll, pos = 1; end

nhyps = length(RESULT.hypnames);
for k = 1:nhyps
    hypname = RESULT.hypnames{k};
    if ~getAll
        pos = strfind(hypname, HYP);
    end
    if ismember(1, pos)
        HV(end+1).name = hypname;
        hypval = RESULT.hypcov(k);
        HV(end).value = hypval;
        
        % Display on screen
        fprintf('%s ==> %g\n', hypname, exp(hypval));
    end
end

end
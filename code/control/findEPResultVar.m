function varIdx = findEPResultVar(sObj, sVar, allObjs, allVars)
%Find a given variable sVar of a given object sObj in the list of result
%variables of E+ read by mlepLoadEPResults(). Cell arrays of all object
%names and all variable names are given by allObjs and allVars.
%
%If sObj is empty then all variables of all objects matching sVar are
%returned.
%
% (C) 2014 by Truong X. Nghiem

if isempty(sObj)
    varIdx = find(strcmpi(sVar, allVars));
    if isempty(varIdx)
        error('Variable %s is not exported in the output file.', sVar);
    end
else
    varIdx = find(strcmpi(sObj, allObjs) & strcmpi(sVar, allVars));
    if isempty(varIdx)
        error('Variable %s of object %s is not exported in the output file.',...
            sVar, sObj);
    end
end

end


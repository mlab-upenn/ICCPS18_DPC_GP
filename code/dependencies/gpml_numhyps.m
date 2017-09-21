function N = gpml_numhyps(F, D)
%GPML_NUMHYPS Returns the number of hyperparameters for a GPML function F.
%   Given a GPML covariance / mean / likelihood function F and the input
%   dimension D, this function returns the total number of hyperparameters.

if ~iscell(F)
    F = {F};
end
N = eval(feval(F{:}));

end


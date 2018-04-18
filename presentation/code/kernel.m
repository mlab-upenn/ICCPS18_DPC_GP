% Define the kernel function
function K = kernel(a, b, param, varargin)
% sqdist = repmat(sum(a.^2,2), [1,size(a,1)]) + repmat(sum(b.^2,2),[1,size(a,1)]) - 2*a*b';
if nargin==4
    K =  varargin{1}*exp(-0.5 * (1/param) * sq_dist(a', b'));
else
    K =  exp(-0.5 * (1/param) * sq_dist(a', b'));
end
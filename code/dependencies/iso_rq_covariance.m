% ISO_RQ_COVARIANCE squared exponential covariance with ARD.
%
% This provides a GPML-compatible covariance function implementing the
% rational quadratic covariance with with isotropic distance measure
% This can be used as a drop-in replacement for covSEiso.
%
% This implementation supports an extended GPML syntax that allows
% calculating the Hessian of K with respect to any pair of
% hyperparameters. The syntax is:
%
%   dK2_didj = iso_rq_covariance(theta, x, z, i, j)
%
% where dK2_didj is \partial^2 K / \partial \theta_i \partial \theta_j,
% and the Hessian is evalauted at K(x, z). As in the derivative API,
% if z is empty, then the Hessian is evaluated at K(x, x).  Note that
% the option of setting z = 'diag' for Hessian computations is not
% supported due to no obvious need.
%
% These Hessians can be used to ultimately compute the Hessian of the
% GP training likelihood (see, for example, exact_inference.m).
%
% The hyperparameters are the same as for covSEiso.
%
% See also COVRQISO, COVFUNCTIONS.

% Copyright (c) 2013--2015 Roman Garnett.

function result = iso_rq_covariance(theta, x, z, i, j)

  % used during gradient and Hessian calculations to avoid constant recomputation
  persistent K;

  % call covRQiso for everything but Hessian calculation
  if (nargin <= 1)
    result = covRQiso;
  elseif (nargin == 2)
    result = covRQiso(theta, x);
  elseif (nargin == 3)
    result = covRQiso(theta, x, z);
  elseif (nargin == 4)
    result = covRQiso(theta, x, z, i);

  % Hessian with respect to \theta_i \theta_j
  else

    % avoid if (isempty(z)) checks
    if (isempty(z))
        z = x;
    end
      
    % ensure i <= j by exploiting symmetry
    if (i > j)
      result = iso_rq_covariance(theta, x, z, j, i);
      return;
    end

    ell = exp(theta(1));
    sf2 = exp(2*theta(2));
    alpha = exp(theta(3));
    D2 = sq_dist(x'/ell,z'/ell);
    
    % Hessian for log(l) and log(l)
    if ((i == 1) && (j == 1))
        K = (1+0.5*D2/alpha);
        result = sf2*(K).^(-alpha-2).*D2.*D2*(1+alpha)/alpha ...
            - 2*sf2*(K).^(-alpha-1).*D2;
    end
    
    % Hessians involving log(sf) and log(ell)
    if (i == 1 && j == 2)
      result = 2 * covRQiso(theta, x, z, i);
      return;
    end
    
    % Hessians involving log(sf) and log(sf)
    if (i == 2 && j == 2)
      result = 2 * covRQiso(theta, x, z, i);
      return;
    end

    % Hessians involving log(sf) and log(alpha)
    if (i == 2 && j == 3)
      result = 2 * covRQiso(theta, x, z, j);
      return;
    end

    % Hessian for log(l) and log(alpha)
    if ((i == 1) && (j == 3))
      result = sf2.*D2.*K.^(-alpha-1).*((1+alpha)/alpha*0.5*D2./K - alpha*log(K));
    end

    % Hessian for log(alpha) and log(alpha)
    if ((i == 3) && (j == 3)) 
       result = sf2*K.^(-alpha).*((0.5*D2./K - alpha*log(K)).^2 + ...
           (0.5*D2).^2/alpha./K^2 - alpha*log(K) + 0.5*D2./K);
    end
    
  end

end
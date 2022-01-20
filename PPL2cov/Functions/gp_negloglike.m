function nlogL = gp_negloglike(xi,sigma,z)

% Negative log-likelihood for the generalized Pareto distribution
% NOTE: Inputs must be same size i.e. size(xi)=size(sigma)=size(z)
% INPUTS:
%   z = threshold exceedances
%   xi = GP shape
%   sigma = GP scale

% some preliminaries
z_norm = z./sigma;
lnsigma = log(sigma);

% calculate likelihood for cases with zero and non-zero shape
nlogL = NaN(size(z));
zeroshp = abs(xi) <= eps;
nlogL(~zeroshp) = lnsigma(~zeroshp) + (1+1./xi(~zeroshp)).*log1p(xi(~zeroshp).*z_norm(~zeroshp));
nlogL(zeroshp) = lnsigma(zeroshp) + z_norm(zeroshp);

% deal with cases where support is violated
bad = xi<0 & z_norm > -1./xi;
nlogL(bad) = Inf;

end
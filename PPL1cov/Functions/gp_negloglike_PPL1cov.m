function pnlogL=gp_negloglike_PPL1cov(param,lambda,const_xi,DATA,GRID)

% Penalised Piecewise-linear GP negative log likelihood
% INPUTS:
%   param = array of GP scale and shape parameters, either
%         = 1x(2*Nnode) array = [sigma(1), ..., sigma(Nnode), xi(1), ..., xi(Nnode)]
%         = 1x(Nnode+1) array = [sigma(1), ..., sigma(Nnode), xi]
%   lambda = 1x2 array of roughness penalties
%          = [lambda_sig, lambda_xi]
%   const_xi = 1x1 logical indicating whether xi is constant
%   DATA = data structure
%   GRID = grid structure
% OUTPUT:
%   pnlogL = penalised negative log likelihood

% Parse inputs
sigma_node=param(1:GRID.Nnode);
xi_node=param(GRID.Nnode+1:end);
if const_xi
    xi_node=0*sigma_node+xi_node;
end
Z=DATA.exceedance.Z;

% check sigma is positive and -0.5 <= xi <= 0
if any(sigma_node<0) || any(xi_node<-0.5) || any(xi_node>0)
    pnlogL=Inf;
    return
end

% Calculate piecewise-linear values of parameters
[sigma,a_sig]=interp_line(DATA.exceedance.X,DATA.exceedance.bin_num,GRID,sigma_node);
[xi,a_xi]=interp_line(DATA.exceedance.X,DATA.exceedance.bin_num,GRID,xi_node);

% check if any data out of range
if any(xi<0 & Z>-sigma./xi)
    pnlogL=Inf;
    return
end

% Negative log-likelihood
nlogL = gp_negloglike(xi,sigma,Z);
if any(abs(imag(nlogL))>0)
    error('imaginary likelihood found')
else
    nlogL=sum(nlogL);
end

% Penalise for roughness (sum of absolute 1st derivatives)
lambda_sig=lambda(1);
pen_sig=sum(abs(a_sig)*lambda_sig);
if const_xi
    pen_xi=0;
else
    lambda_xi=lambda(2);
    pen_xi=sum(abs(a_xi)*lambda_xi);
end

pnlogL = nlogL + pen_xi + pen_sig;

end
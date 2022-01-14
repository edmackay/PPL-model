function DATA=set_initial_GP_voronoi_1cov(DATA)

% Assumes input structure with fields:
%   DATA.exceedance.X - Nx1 array of covariates
%   DATA.exceedance.Y - Nx1 array of threshold exceedances
%   DATA.exceedance.bin_num - Nx1 array of bin_numbers for each exceedance
%   DATA.exceedance.nearest_node - Nx1 array of nearest nodes for each exceedance

Nnode=max(DATA.exceedance.bin_num);
param_indep_bin=NaN(Nnode,2);

% calculate independent ML estimates in each bin
for i = 1:Nnode
    % set bins
    bin=DATA.exceedance.nearest_node==i;
    z=DATA.exceedance.Z(bin);
    
    % calculate ML estimators of GP params
    [sigma, xi]=gp_fit_constrained(z,-0.5,0);
    param_indep_bin(i,:)=[sigma,xi];
end
sig1=param_indep_bin(:,1);
xi1=param_indep_bin(:,2);

% set initial guess for sigma as independent moment estimates and initial guess for xi as just below 0
xi0=-0.01;
param0=[param_indep_bin(:,1); xi0];

% optimise piecewise-constant estimates
param_const_xi=fminsearch(@(p)gplike_PC_constxi(p,DATA.exceedance.Z,DATA.exceedance.nearest_node),param0);
sig0=param_const_xi(1:Nnode);
xi0=param_const_xi(end);

% set outputs
DATA.voronoi.xi_const.sigma=sig0;
DATA.voronoi.xi_const.xi=xi0;
DATA.voronoi.xi_vary.sigma=sig1;
DATA.voronoi.xi_vary.xi=xi1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function nlogL=gplike_PC_constxi(param,Z,A)
        Sig=param(1:end-1);
        Xi=param(end);
        
        % check sigma is positive and -0.5 <= xi <= 0
        if any(Sig<0) || Xi<-0.5 || Xi>0
            nlogL=Inf;
            return
        end
        
        % lookup sigma and xi for each value of z
        xi_dat=0*Z+Xi;
        sig_dat=Sig(A);
        
        % calculate negative log-likelihood
        nlogL_dat = gp_negloglike(xi_dat,sig_dat,Z);
        nlogL = sum(nlogL_dat);
    end

end
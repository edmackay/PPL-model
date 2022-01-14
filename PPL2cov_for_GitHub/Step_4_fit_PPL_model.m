clear

addpath('Functions')

filein='NorthSea_crossval_irregular_4_pensig_1_penxi_0.mat';

load(filein,'DATA','GRID','CV')

P=0.99;                             % conditional quantiles
Npoints=10*DATA.exceedance.Ndat;    % Number of points to simulate
xrange=[2 18];                      % x range for exceedance plots
yrange=[1e-4 1];                    % y range for exceedance plots
Nboot=100;                          % number of bootstrap trials

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       EDIT INPUTS ONLY ABOVE THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fit model
param_nodes=fit_PPL_model_2cov(DATA,GRID,CV,Nboot);

% save data
fileout=strrep(filein,'crossval','fitted');
save(fileout,'DATA','GRID','CV','param_nodes')

% parameter estimates
if CV.setup.penshape==0
    const_xi=true;
else
    const_xi=false;
end
plot_parameter_estimates_2cov(DATA,GRID,param_nodes,const_xi)

% conditional quantiles
plot_PPL_quantile_2cov(DATA,GRID,param_nodes,P)

% simulate under model
DATA_mod=simulate_model_2cov(DATA,GRID,param_nodes,Npoints);

% exceedance plots by sector 
plot_exceedance_by_sector_2cov(DATA,DATA_mod,xrange,yrange)

% exceedance plots by bin
plot_exceedance_by_bin_2cov(DATA,DATA_mod,GRID,xrange,yrange)
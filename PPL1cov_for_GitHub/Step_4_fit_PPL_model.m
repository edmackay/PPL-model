clear

addpath('Functions')

filein='NorthSea_Season_crossval_4node_penxi_0.mat'; type='season';

load(filein,'DATA','GRID','CV')

P=[0.9 0.99 0.999];                 % conditional quantiles
Npoints=100*DATA.exceedance.Ndat;   % Number of points to simulate
Nboot=100;                          % number of bootstrap trials
xrange=[2 18];                      % x range for exceedance plots
yrange=[1e-4 1];                    % y range for exceedance plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       EDIT INPUTS ONLY ABOVE THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fit model
param_nodes=fit_PPL_model_1cov(DATA,GRID,CV,Nboot);

% save data
fileout=strrep(filein,'crossval','fitted');
save(fileout,'DATA','GRID','CV','param_nodes')

% parameter estimates
if CV.setup.penshape==0
    const_xi=true;
else
    const_xi=false;
end
plot_parameter_estimates_1cov(DATA,GRID,param_nodes,const_xi)

% conditional quantiles
figure
plot_PPL_quantile_1cov(DATA,GRID,param_nodes,P)

% simulate under model
DATA_mod=simulate_model_1cov(DATA,GRID,param_nodes,Npoints);

% exceedance plots by sector 
figure
plot_exceedance_by_sector_1cov(DATA,DATA_mod,type,xrange,yrange)

% exceedance plots by bin
plot_exceedance_by_bin_1cov(DATA,DATA_mod,GRID,xrange,yrange)
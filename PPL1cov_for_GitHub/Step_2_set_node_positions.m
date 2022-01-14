clear

addpath('Functions')

% Load data
load('NorthSea_Season_initialised.mat')

% seasonal node locations
nodes=20:90:360;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       EDIT INPUTS ONLY ABOVE THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save output file
[DATA, GRID]=create_grid(DATA,nodes);
fileout=[DATA.name.dataset '_binned_' num2str(length(nodes)) 'node'];

% calculate initial GP parameter estimates from Voronoi partition
DATA=set_initial_GP_voronoi_1cov(DATA);

% check grid
plot_grid_checks_1cov(DATA,GRID)

save(fileout,'DATA','GRID')



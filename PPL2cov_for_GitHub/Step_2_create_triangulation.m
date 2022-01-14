clear

addpath('Functions')
set(0,'defaultfigurewindowstyle','docked')

% Load data
load('NorthSea_initialised.mat')

% % create regular grid
% regular_grid=true;              % create a regular grid
% nodes.x1=[120 210 270];         % list of marginal node positions in X1
% nodes.x2=[20 200];              % list of marginal node positions in X2

% create irregular grid
regular_grid=false;             % create an irregular grid
nodes=[270 40                   % N x 2 list of node positions in X1 and X2
       165 240
       270 200
       50 180];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       EDIT INPUTS ONLY ABOVE THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create grid
[DATA, GRID, fileout]=create_grid(DATA,regular_grid,nodes);

% calculate initial GP parameter estimates from Voronoi partition
DATA=set_initial_GP_voronoi_2cov(DATA);

% save output file
save(fileout,'DATA','GRID')

% Check triangulation
plot_grid_checks_2cov(DATA,GRID)
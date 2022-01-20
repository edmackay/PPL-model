clear

addpath('Functions')
set(0,'defaultfigurewindowstyle','docked')

% Load data
load('North_Sea_initialised.mat')

% % create regular grid
% regular_grid=true;              % create a regular grid
% nodes.x1=[120 210 270];         % list of marginal node positions in X1
% nodes.x2=[20 200];              % list of marginal node positions in X2

% create irregular grid
regular_grid=false;                       % create an irregular grid
nodes.x1=[270, 165, 170, 270,  50,  50];  % list of node positions in X1
nodes.x2=[ 40, 240,  80, 200, 180, 300];  % list of node positions in X2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       EDIT INPUTS ONLY ABOVE THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create grid
[DATA, GRID, fileout]=create_grid(DATA,regular_grid,nodes);

% calculate initial GP parameter estimates from Voronoi partition
DATA=set_initial_GP_voronoi_2cov(DATA);

% Check triangulation
plot_grid_checks_2cov(DATA,GRID)

% save output file
save(fileout,'DATA','GRID')
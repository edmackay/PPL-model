clear

addpath('Functions')

% Load data
load('North_Sea_Data.mat')

% set input data structure
DATA.X=Direction;                    % Ndat x 1 array of covariate values
DATA.Y=Hs_peak;                      % Ndat x 1 array of response values
DATA.name.dataset='North_Sea_directional';   % Name of dataset
DATA.name.X='Direction';             % Name of covariate
DATA.name.Y='Sig. wave height';      % Name of response
DATA.unit.X='deg';                   % Unit of covariate
DATA.unit.Y='m';                     % Unit of response
DATA.Nyears=54;                      % Length of dataset in years

% settings for gridded estimates
DATA.grid.x=0:1:360;                 % Resolution of grid
DATA.grid.Nnearest=50;               % Number of nearest points for threshold and GP parameter estimates
DATA.grid.KD_bw=5;                   % Bandwidth for kernel density estimate of covariate PDF
DATA.grid.Thresh_bw=5;               % Bandwidth for smoothing threshold estimate
DATA.thresh_NEP=0.8;                 % Threshold non-exceedance probability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       EDIT INPUTS ONLY ABOVE THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate initial gridded estimates
DATA=set_initial_gridded_estimates_1cov(DATA);

% visualise gridded estimates
figure
subplot(2,2,1)
plot_initial_gridded_estimates_1cov(DATA,'pdf')

subplot(2,2,2)
plot_initial_gridded_estimates_1cov(DATA,'thresh')

subplot(2,2,3)
plot_initial_gridded_estimates_1cov(DATA,'sigma')

subplot(2,2,4)
plot_initial_gridded_estimates_1cov(DATA,'xi')

% save initialised data
save([DATA.name.dataset '_initialised'],'DATA')


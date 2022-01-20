clear

addpath('Functions')

% Load data
load('North_Sea_Data.mat')

% set input data structure
DATA.X=[Direction, Season];       % Ndat x 2 array of covariate values
DATA.Y=Hs_peak;                   % Ndat x 1 array of response values
DATA.name.dataset='North_Sea11';    % Name of dataset
DATA.name.X1='Direction';         % Name of 1st covariate
DATA.name.X2='Season';            % Name of 2nd covariate
DATA.name.Y='Sig. wave height';   % Name of response
DATA.unit.X1='deg';               % Unit of 1st covariate
DATA.unit.X2='day';               % Unit of 2nd covariate
DATA.unit.Y='m';                  % Unit of response
DATA.Nyears=54;                   % Length of dataset in years

% settings for gridded estimates
DATA.grid.x1=0:2:360;             % Resolution of grid for 1st covariate
DATA.grid.x2=0:2:360;             % Resolution of grid for 2nd covariate
DATA.grid.Nnearest=50;            % Number of nearest points for threshold and initial GP parameter estimates
DATA.grid.KD_bw=10;               % Bandwidth for kernel density estimate of covariate joint PDF
DATA.grid.Thresh_bw=10;           % Bandwidth for Gaussian smoothing kernel for threshold
DATA.grid.sigma_bw=10;            % Bandwidth for Gaussian smoothing kernel for scale
DATA.thresh_NEP=0.7;              % Threshold non-exceedance probability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       EDIT INPUTS ONLY ABOVE THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate initial gridded estimates
DATA=set_initial_gridded_estimates_2cov(DATA);

% visualise gridded estimates
figure
subplot(2,2,1)
plot_initial_gridded_estimates_2cov(DATA,'pdf',1)

subplot(2,2,2)
plot_initial_gridded_estimates_2cov(DATA,'thresh',1)

subplot(2,2,3)
plot_initial_gridded_estimates_2cov(DATA,'sigma',1)

subplot(2,2,4)
plot_initial_gridded_estimates_2cov(DATA,'xi',1)

% save initialised data
save([DATA.name.dataset '_initialised'],'DATA')


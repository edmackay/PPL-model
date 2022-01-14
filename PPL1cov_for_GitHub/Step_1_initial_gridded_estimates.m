clear

addpath('Functions')

% Load data
load('NorthSeaData.mat','Pks')

% set input data structure
DATA.X=Pks.X(:,2);                   % Ndat x 2 array of covariate values
DATA.Y=Pks.Y;                        % Ndat x 1 array of response values
DATA.name.dataset='NorthSea_Season'; % Name of dataset
DATA.name.X='Season';                % Name of 1st covariate
DATA.name.Y='Sig. wave height';      % Name of response
DATA.unit.X='day';                   % Unit of 1st covariate
DATA.unit.Y='m';                     % Unit of response
DATA.Nyears=60;                      % Length of dataset in years

% settings for gridded estimates
DATA.grid.x=0:1:360;                 % Resolution of grid for 1st covariate
DATA.grid.Nnearest=100;              % Number of nearest points for threshold and GP estimates
DATA.grid.Thresh_binwidth=15;        % Binwidth for smoothing threshold estimate
DATA.grid.KD_bw=15;                  % Bandwidth for kernel density estimate of covariate joint PDF
DATA.thresh_NEP=0.7;                 % Threshold non-exceedance probability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       EDIT INPUTS ONLY ABOVE THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate initial gridded estimates
DATA=set_initial_gridded_estimates_1cov(DATA);

% visualise gridded estimates
figure
plot_initial_gridded_estimates_1cov(DATA,'pdf')

figure
plot_initial_gridded_estimates_1cov(DATA,'thresh')

figure
plot_initial_gridded_estimates_1cov(DATA,'sigma')

figure
plot_initial_gridded_estimates_1cov(DATA,'xi')

% save initialised data
save([DATA.name.dataset '_initialised'],'DATA')


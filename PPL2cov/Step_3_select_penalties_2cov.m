clear

addpath('Functions')

filein='North_Sea_binned_irregular_6.mat';

load(filein,'GRID','DATA')

% Cross validation setup
CV.setup.Nrepeats=5;    % number of repeats
CV.setup.Ngroups=5;     % number of cross-validation groups
CV.setup.Npenalty=10;   % number of values tried for each lambda
CV.setup.penLB=-1;      % lower bound (log10) for lambda range
CV.setup.penUB=5;       % upper bound (log10) for lambda range
CV.setup.penshape=0;    % 0 - constant shape (i.e., inifinite roughness penalty), 1 - one roughness penalty on shape
CV.setup.penscale=1;    % 1 - one roughness penalty on scale, 2 - a roughness penalty on scale for each covariate
CV.options=optimoptions(@fmincon, 'Display','notify','MaxFunctionEvaluations',1e5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       EDIT INPUTS ONLY ABOVE THIS LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run cross-validation
CV = crossval_2cov(DATA,GRID,CV);

plot_cross_validation_results_2cov(CV)

% save cross validation data
fileout=strrep(filein,'binned','crossval');
fileout=[fileout(1:end-4) '_pensig_' num2str(CV.setup.penscale) '_penxi_' num2str(CV.setup.penshape) '.mat'];
save(fileout,'DATA','GRID','CV')


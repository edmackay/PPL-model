function CV=crossval_2cov(DATA,GRID,CV)

% parse inputs
Ndat=DATA.exceedance.Ndat;
Nnode=GRID.Nnode_unique;
Npen=CV.setup.Npenalty;
Nrep=CV.setup.Nrepeats;
Ngrp=CV.setup.Ngroups;
penshape=CV.setup.penshape;
penscale=CV.setup.penscale;

% Create grid of lambda values to search
lambdavals=logspace(CV.setup.penLB,CV.setup.penUB,Npen);
if penscale == 1 && penshape == 0
    % Case A: common penalty for sigma, constant xi
    lambda = [lambdavals(:), lambdavals(:)];
    Nlambda=Npen;
    const_xi=true;
    param0=[DATA.voronoi.xi_const.sigma; DATA.voronoi.xi_const.xi];
elseif penscale == 2 && penshape==0
    % Case B: different penalties for sigma, constant xi
    [lambda_sig1, lambda_sig2] = meshgrid(lambdavals);
    lambda = [lambda_sig1(:), lambda_sig2(:)];
    Nlambda=Npen^2;
    const_xi=true;
    param0=[DATA.voronoi.xi_const.sigma; DATA.voronoi.xi_const.xi];
elseif penscale == 1 && penshape == 1
    % Case C: common penalty for sigma, common penalty for xi
    [lambda_sig, lambda_xi] = meshgrid(lambdavals);
    lambda = [lambda_sig(:), lambda_sig(:), lambda_xi(:), lambda_xi(:)];
    Nlambda=Npen^2;
    const_xi=false;
    param0=[DATA.voronoi.xi_vary.sigma; DATA.voronoi.xi_vary.xi];
elseif penscale == 2 && penshape == 1
    % Case D: different penalties for sigma, common penalty for xi
    [lambda_sig1, lambda_sig2, lambda_xi] = meshgrid(lambdavals);
    lambda = [lambda_sig1(:), lambda_sig2(:), lambda_xi(:), lambda_xi(:)];
    Nlambda=Npen^3;
    const_xi=false;
    param0=[DATA.voronoi.xi_vary.sigma; DATA.voronoi.xi_vary.xi];
else 
    error('Only cases A-D supported')
end

% set any zero shape to slightly negative
for i=Nnode+1:length(param0)
    if param0(i)==0
        param0(i)=-1e-6;
    end
end

% Partition exceedances into groups for each repetition
CV.groups = zeros(Ndat,Nrep);
for i=1:Nrep
    CV.groups(:,i) = CVgroups(Ngrp,DATA.exceedance.bin_num_unique);
end

% Set upper and lower bounds for fmincon (sigma>0 and -0.5<xi<0)
param_LB=NaN*param0;
param_LB(1:Nnode)=0;
param_LB(Nnode+1:end)=-0.5;
param_UB=NaN*param0;
param_UB(1:Nnode)=10*max(DATA.exceedance.Z);
param_UB(Nnode+1:end)=0;

% Cross Validation procedure
predictive_likelihood=zeros(Nlambda,Nrep);
for irep=1:Nrep
    for igrp=1:Ngrp
        display(['   Repeat ' num2str(irep) ', Group ' num2str(igrp)])
        
        % set samples to fit and to validate
        ind_val = CV.groups(:,irep)==igrp;
        
        DATA_val.exceedance.X=DATA.exceedance.X(ind_val,:);
        DATA_val.exceedance.Z=DATA.exceedance.Z(ind_val,:);
        DATA_val.exceedance.bin_num=DATA.exceedance.bin_num(ind_val,:);
        DATA_val.exceedance.Ndat=sum(ind_val);
        
        DATA_fit.exceedance.X=DATA.exceedance.X(~ind_val,:);
        DATA_fit.exceedance.Z=DATA.exceedance.Z(~ind_val,:);
        DATA_fit.exceedance.bin_num=DATA.exceedance.bin_num(~ind_val,:);
        DATA_fit.exceedance.Ndat=sum(~ind_val);
        
        % Constrained optimisation
        fprintf('         Penalty');
        for ilambda=1:Nlambda
            fprintf(' %d', ilambda);
            if mod(ilambda,25)==0
                fprintf('\n                ');
            end
            
            % check initial log-likelihood
            init=gp_negloglike_PPL2cov(param0,[0,0,0,0],const_xi,DATA_fit,GRID);
            if isinf(init)
                while isinf(init)
                    param0(1:Nnode)=1.5*param0(1:Nnode);
                    init=gp_negloglike_PPL2cov(param0,[0,0,0,0],const_xi,DATA_fit,GRID);
                end
            end
            
            % optimise parameters for included data
            paramhat=fmincon(@(p)gp_negloglike_PPL2cov(p,lambda(ilambda,:),const_xi,DATA_fit,GRID),param0,[],[],[],[],param_LB,param_UB,[],CV.options);

            % Compute likelihood for excluded data (with zero penalty)
            val_like=gp_negloglike_PPL2cov(paramhat,[0,0,0,0],const_xi,DATA_val,GRID);
            
            % update predictive likelihood
            predictive_likelihood(ilambda,irep) = predictive_likelihood(ilambda,irep) + val_like;
        end
        fprintf('\n');
    end
end

% Take 'optimal' parameter values to be the combination leading to the smallest mean cross validation score across all repeats.
mean_pred_like= mean(predictive_likelihood,2);
[CV.results.predictive_likelihood_min, optlambdaind] = min(mean_pred_like);
CV.results.predictive_likelihood_min_jackknife_range=max(jackknife(@mean,predictive_likelihood(optlambdaind,:)))-min(jackknife(@mean,predictive_likelihood(optlambdaind,:)));
CV.results.lambda_optimal = lambda(optlambdaind,:);

% Rough rule to decide on acceptable lambda combinations is allowing the mean predictive performance to be
CV.results.predictive_likelihood_accepted = CV.results.predictive_likelihood_min + CV.results.predictive_likelihood_min_jackknife_range;
CV.results.lambda_accepted= lambda(mean_pred_like<CV.results.predictive_likelihood_accepted,:);

% reshape predictive likelihood so it is easier to interpret
% PL_out = Nrep x Nlambda1 x Nlambda2 x ...
if penscale == 1 && penshape == 0
    % Case A: common penalty for sigma, constant xi
    PL_out = predictive_likelihood';
elseif penscale == 2 && penshape==0
    % Case B: different penalties for sigma, constant xi
    PL_out=zeros(Nrep,Npen,Npen);
    for i=1:Nrep
        PL_out(i,:,:) = reshape(predictive_likelihood(:,i),Npen,Npen);
    end
elseif penscale == 1 && penshape == 1
    % Case C: common penalty for sigma, common penalty for xi
    PL_out=zeros(Nrep,Npen,Npen);
    for i=1:Nrep
        PL_out(i,:,:) = reshape(predictive_likelihood(:,i),Npen,Npen);
    end
elseif penscale == 2 && penshape == 1
    % Case D: different penalties for sigma, common penalty for xi
    PL_out=zeros(Nrep,Npen,Npen,Npen);
    for i=1:Nrep
        PL_out(i,:,:,:) = reshape(predictive_likelihood(:,i),Npen,Npen,Npen);
    end
end
CV.results.predictive_likelihood = PL_out;
CV.results.lambda=lambdavals;
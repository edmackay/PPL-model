function CV=crossval_1cov(DATA,GRID,CV)

% parse inputs
Ndat=DATA.exceedance.Ndat;
Nnode=GRID.Nnode;
Npen=CV.setup.Npenalty;
Nrep=CV.setup.Nrepeats;
Ngrp=CV.setup.Ngroups;
penshape=CV.setup.penshape;

% Create grid of lambda values to search
lambdavals=logspace(CV.setup.penLB,CV.setup.penUB,Npen);
if penshape == 0
    % Case A: constant shape
    lambda = lambdavals(:);
    Nlambda=Npen;
    const_xi=true;
    param0=[DATA.voronoi.xi_const.sigma; DATA.voronoi.xi_const.xi];
elseif penshape == 1
    % Case B: variable shape
    [lambda_sig, lambda_xi] = meshgrid(lambdavals);
    lambda = [lambda_sig(:), lambda_xi(:)];
    Nlambda=Npen^2;
    const_xi=false;
    xi0=DATA.voronoi.xi_vary.xi;
    xi0(xi0==0)=-0.001;
    param0=[DATA.voronoi.xi_vary.sigma; xi0];
end

% Partition exceedances into groups for each repetition
CV.groups = zeros(Ndat,Nrep);
for i=1:Nrep
    CV.groups(:,i) = CVgroups(Ngrp,DATA.exceedance.bin_num);
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
            
            % optimise parameters for included data
            paramhat=fmincon(@(p)gp_negloglike_PPL1cov(p,lambda(ilambda,:),const_xi,DATA_fit,GRID),param0,[],[],[],[],param_LB,param_UB,[],CV.options);
            
            % Compute likelihood for excluded data (with zero penalty)
            val_like=gp_negloglike_PPL1cov(paramhat,[0,0],const_xi,DATA_val,GRID);
            
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
if penshape == 0
    % Case A: constant xi
    PL_out = predictive_likelihood';
elseif penshape == 1
    % Case B: variable xi
    PL_out=zeros(Nrep,Npen,Npen);
    for i=1:Nrep
        PL_out(i,:,:) = reshape(predictive_likelihood(:,i),Npen,Npen);
    end
end
CV.results.predictive_likelihood = PL_out;
CV.results.lambda=lambdavals;
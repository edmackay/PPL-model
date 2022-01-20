function param_nodes=fit_PPL_model_1cov(DATA,GRID,CV,Nboot)

lambda_opt=CV.results.lambda_optimal;
% lambda_opt=0.1;
if CV.setup.penshape == 0
    const_xi=true;
    param0=[DATA.voronoi.xi_const.sigma; DATA.voronoi.xi_const.xi];
elseif CV.setup.penshape == 1
    const_xi=false;
    xi0=DATA.voronoi.xi_vary.xi;
    xi0(xi0==0)=-0.001;
    param0=[DATA.voronoi.xi_vary.sigma; xi0];
end

% Set upper and lower bounds for fmincon (sigma>0 and -0.5<xi<0)
Nnode=GRID.Nnode;
param_LB=NaN*param0;
param_LB(1:Nnode)=0;
param_LB(Nnode+1:end)=-0.5;
param_UB=NaN*param0;
param_UB(1:Nnode)=10*max(DATA.exceedance.Z);
param_UB(Nnode+1:end)=0;

param_nodes=zeros(Nnode,2,Nboot);
if Nboot==1
    % fit model
    paramhat=fmincon(@(p)gp_negloglike_PPL1cov(p,lambda_opt,const_xi,DATA,GRID),param0,[],[],[],[],param_LB,param_UB,[],CV.options);
    param_nodes(:,1)=paramhat(1:Nnode);
    param_nodes(:,2)=paramhat(Nnode+1:end);
else
    Ndat=DATA.exceedance.Ndat;
    for i=1:Nboot
        fprintf('Fitting model for bootstrap trial %d\n',i)
        % resample data
        inds=ceil(Ndat*rand(Ndat,1));
        DATAboot.exceedance.X=DATA.exceedance.X(inds);
        DATAboot.exceedance.Z=DATA.exceedance.Z(inds);
        DATAboot.exceedance.bin_num=DATA.exceedance.bin_num(inds);
        paramhat=fmincon(@(p)gp_negloglike_PPL1cov(p,lambda_opt,const_xi,DATAboot,GRID),param0,[],[],[],[],param_LB,param_UB,[],CV.options);
        param_nodes(:,1,i)=paramhat(1:Nnode);
        param_nodes(:,2,i)=paramhat(Nnode+1:end);
    end
end
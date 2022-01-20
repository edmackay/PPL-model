function DATA_mod=simulate_model_2cov(DATA,GRID,param_nodes,Npoints)

% param_nodes = Nnode x 2 list of sigma and xi at unique node locations

x1=DATA.grid.x1;
x2=DATA.grid.x2;
dx1=x1(2)-x1(1);
dx2=x2(2)-x2(1);
PDF=DATA.grid.pdf;

% normalise data
PDF=PDF/sum(PDF(:)*dx1*dx2);

% marginal distribution of X1
Fx1=cumsum(sum(PDF,1))*dx1*dx2;
Fx1(1)=0;
Fx1(end)=1;

% conditional distribution of x2 given x1
Fx2gx1=cumsum(PDF,1);
Fx2gx1=Fx2gx1./repmat(Fx2gx1(end,:),length(x2),1);
Fx2gx1(1,:)=0;

% find conditional inverse
nv=500;
vv=linspace(0,1,nv);
Finv=zeros(nv,length(x1));
for i=1:length(x1)
    Finv(:,i)=interp1(Fx2gx1(:,i),x2,vv);
end

% simulate for each bootstrap
Nboot=size(param_nodes,3);

DATA_mod.X=zeros(Npoints,2,Nboot);
DATA_mod.Y=zeros(Npoints,Nboot);
DATA_mod.Z=zeros(Npoints,Nboot);
DATA_mod.U=zeros(Npoints,Nboot);
DATA_mod.bin_num=zeros(Npoints,Nboot);

for j=1:Nboot
    if Nboot>1
        disp(['Simulating for bootstrap trial ' num2str(j)])
    end
    
    % simulate covariate data
    u=rand(Npoints,1);
    v=rand(Npoints,1);
    X1=interp1(Fx1,x1,u);
    X2=interp2(x1,vv,Finv,X1,v);
    U=interp2(x1,x2,DATA.grid.thresh,X1,X2);
    
    % shift data if grid is regular
    if GRID.isregular
        shiftx1 = X1<min(GRID.nodes.position(:,1));
        shiftx2 = X2<min(GRID.nodes.position(:,2));
        X1(shiftx1)=X1(shiftx1)+360;
        X2(shiftx2)=X2(shiftx2)+360;
    end
    X=[X1,X2];
    
    % calculate bin numbers
    bin_num=zeros(Npoints,1);
    for i=1:GRID.Ntri_repeated
        in = inpolygon(X1,X2,GRID.nodes.position_repeated(GRID.tri.indices(i,:),1),GRID.nodes.position_repeated(GRID.tri.indices(i,:),2));
        bin_num(in)=i;
    end
    
    % interpolate GP parameters
    SIGMA=interp_tri(X,bin_num,GRID,param_nodes(:,1,j));
    XI=interp_tri(X,bin_num,GRID,param_nodes(:,2,j));
    
    % simulate GP random variables
    Y=gprnd(XI,SIGMA,U);

    % set outputs
    DATA_mod.X(:,:,j)=X;
    DATA_mod.Y(:,j)=Y;
    DATA_mod.Z(:,j)=Y-U;
    DATA_mod.U(:,j)=U;
    DATA_mod.bin_num(:,j)=bin_num;
end

end
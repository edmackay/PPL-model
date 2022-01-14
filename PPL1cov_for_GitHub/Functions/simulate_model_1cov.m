function DATA_mod=simulate_model_1cov(DATA,GRID,param_nodes,Npoints)

% param_nodes = Nnode x 2 list of sigma and xi at unique node locations

Nboot=size(param_nodes,3);

DATA_mod.X=zeros(Npoints,Nboot);
DATA_mod.Y=zeros(Npoints,Nboot);
DATA_mod.Z=zeros(Npoints,Nboot);
DATA_mod.U=zeros(Npoints,Nboot);
DATA_mod.bin_num=zeros(Npoints,Nboot);

for i=1:Nboot
    if Nboot>1
        disp(['Simulating for bootstrap trial ' num2str(i)])
    end
    
    x=DATA.grid.x;
    dx=x(2)-x(1);
    PDF=DATA.grid.pdf;
    
    % normalise data
    PDF=PDF/sum(PDF*dx);
    
    % marginal distribution of X1
    Fx1=cumsum(PDF)*dx;
    Fx1(1)=0;
    Fx1(end)=1;
    
    % simulate covariate data
    p=rand(Npoints,1);
    X=interp1(Fx1,x,p);
    U=interp1(x,DATA.grid.thresh,X);
    
    % shift data if grid is regular
    shift = X<min(GRID.nodes.position(:,1));
    X(shift)=X(shift)+360;
    
    % calculate bin numbers
    bin_num=zeros(Npoints,1);
    nodes=[GRID.nodes.position; GRID.nodes.position(1)+360];
    for j=1:GRID.Nnode
        in = X>=nodes(j) & X<=nodes(j+1);
        bin_num(in)=j;
    end
    
    % interpolate GP parameters
    SIGMA=interp_line(X,bin_num,GRID,param_nodes(:,1,i));
    XI=interp_line(X,bin_num,GRID,param_nodes(:,2,i));
    
    % simulate GP random variables
    Y=gprnd(XI,SIGMA,U);
    Z=Y-U;
    
    % set outputs
    DATA_mod.X(:,i)=X;
    DATA_mod.Y(:,i)=Y;
    DATA_mod.Z(:,i)=Z;
    DATA_mod.U(:,i)=U;
    DATA_mod.bin_num(:,i)=bin_num;
end
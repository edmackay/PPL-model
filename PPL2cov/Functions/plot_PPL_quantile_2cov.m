function plot_PPL_quantile_2cov(DATA,GRID,param_nodes,P)

% get threshold and grid
u=DATA.grid.thresh;
x1=DATA.grid.x1;
x2=DATA.grid.x2;
[X1,X2]=meshgrid(x1,x2);

% shift grid if regular
if GRID.isregular
    X1(X1<min(GRID.nodes.position(:,1)))=X1(X1<min(GRID.nodes.position(:,1)))+360;
    X2(X2<min(GRID.nodes.position(:,2)))=X2(X2<min(GRID.nodes.position(:,2)))+360;
end

% set bin numbers
bin_num_x=zeros(size(X1));
for i=1:GRID.Ntri_repeated
    in = inpolygon(X1(:),X2(:),GRID.nodes.position_repeated(GRID.tri.indices(i,:),1),GRID.nodes.position_repeated(GRID.tri.indices(i,:),2));
    bin_num_x(in)=i;
end

% set parameter values
if size(param_nodes,3)>1
    param_nodes=mean(param_nodes,3);
end
sigma=interp_tri([X1(:), X2(:)],bin_num_x(:),GRID,param_nodes(:,1));
xi=interp_tri([X1(:), X2(:)],bin_num_x(:),GRID,param_nodes(:,2));
sigma=reshape(sigma,length(x2),length(x1));
xi=reshape(xi,length(x2),length(x1));

% calculate quantiles
Q=gpinv(P,xi,sigma,u);

figure
hold on; grid on
surf(x1,x2,Q,Q,'EdgeColor','none')
% plot3(DATA.X(:,1),DATA.X(:,2),DATA.Y,'o','color',[1 1 1]*0.6)
xlim([0 360])
xticks(0:90:360)
ylim([0 360])
yticks(0:90:360)
zlabel([DATA.name.Y ' [' DATA.unit.Y ']'])
xlabel([DATA.name.X1 ' [' DATA.unit.X1 ']'])
ylabel([DATA.name.X2 ' [' DATA.unit.X2 ']'])
light('Position',[0 0 3*max(DATA.Y)],'Style','local')
lighting gouraud
view(3)
rotate3d on
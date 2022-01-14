function plot_parameter_estimates_2cov(DATA,GRID,param_nodes,const_xi)

set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')

Nboot=size(param_nodes,3);

sigma=squeeze(param_nodes(:,1,:));
xi=squeeze(param_nodes(:,2,:));

sigmean=mean(sigma,2);
ximean=mean(xi,2);

sigstd=std(sigma,[],2);
xistd=std(xi,[],2);

figure
trisurf(GRID.tri.indices,GRID.nodes.position_repeated(:,1),GRID.nodes.position_repeated(:,2),sigmean(GRID.nodes.number),'facecolor','interp')
xlabel([DATA.name.X1 ' [' DATA.unit.X1 ']'])
ylabel([DATA.name.X2 ' [' DATA.unit.X2 ']'])
zlabel('mean$(\sigma(x))$')

if Nboot>1
    figure
    trisurf(GRID.tri.indices,GRID.nodes.position_repeated(:,1),GRID.nodes.position_repeated(:,2),sigstd(GRID.nodes.number),'facecolor','interp')
    xlabel([DATA.name.X1 ' [' DATA.unit.X1 ']'])
    ylabel([DATA.name.X2 ' [' DATA.unit.X2 ']'])
    zlabel('STD$(\sigma(x))$')
end

if const_xi && Nboot>1
    figure
    h=histogram(xi(1,:),'Normalization','pdf');
    h.FaceColor=[1 1 1]*0.5;
    xlabel('GP shape')
    ylabel('Probability density')
elseif ~const_xi
    figure
    trisurf(GRID.tri.indices,GRID.nodes.position_repeated(:,1),GRID.nodes.position_repeated(:,2),ximean(GRID.nodes.number),'facecolor','interp')
    xlabel([DATA.name.X1 ' [' DATA.unit.X1 ']'])
    ylabel([DATA.name.X2 ' [' DATA.unit.X2 ']'])
    zlabel('mean$(\xi(x))$')
    
    figure
    trisurf(GRID.tri.indices,GRID.nodes.position_repeated(:,1),GRID.nodes.position_repeated(:,2),xistd(GRID.nodes.number),'facecolor','interp')
    xlabel([DATA.name.X1 ' [' DATA.unit.X1 ']'])
    ylabel([DATA.name.X2 ' [' DATA.unit.X2 ']'])
    zlabel('STD$(\xi(x))$')
end
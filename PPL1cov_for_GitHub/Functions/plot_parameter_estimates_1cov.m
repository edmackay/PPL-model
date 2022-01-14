function plot_parameter_estimates_1cov(DATA,GRID,param_nodes,const_xi)

nodes=GRID.nodes.position([end-1, 1:end]);
nodes(1)=nodes(1)-360;
Nboot=size(param_nodes,3);

sigma=squeeze(param_nodes([end 1:end 1],1,:));
xi=squeeze(param_nodes([end 1:end 1],2,:));

figure
hold on; box on
plot(nodes,sigma,'color',[1 1 1]*0.6)
plot(nodes,mean(sigma,2),'k-o','LineWidth',2)
xlabel([DATA.name.X ' [' DATA.unit.X ']'])
ylabel('GP scale')
xlim([0 360])
xticks(0:45:360)

if const_xi && Nboot>1
    figure
    hold on; box on
    h=histogram(squeeze(param_nodes(1,2,:)),'Normalization','pdf');
    h.FaceColor=[1 1 1]*0.6;
    xlabel('GP shape')
    ylabel('Probability density')
elseif ~const_xi
    figure
    hold on; box on
    plot(nodes,xi,'color',[1 1 1]*0.6)
    plot(nodes,mean(xi,2),'k-o','LineWidth',2)
    xlabel([DATA.name.X ' [' DATA.unit.X ']'])
    ylabel('GP shape')
    xlim([0 360])
    xticks(0:45:360)
end
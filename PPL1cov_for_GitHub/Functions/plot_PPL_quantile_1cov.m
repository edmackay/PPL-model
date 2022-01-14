function plot_PPL_quantile_1cov(DATA,GRID,param_nodes,P)

if size(param_nodes,3)>1
    param_nodes=mean(param_nodes,3);
end

% interpolate parameters to regular intervals
x=DATA.grid.x(:);
u=DATA.grid.thresh(:);
bin_num_x=zeros(length(x),1);
x(x<GRID.nodes.position(1))=x(x<GRID.nodes.position(1))+360;
nodes=[GRID.nodes.position; GRID.nodes.position(1)+360];
for i=1:GRID.Nnode
    in = x>=nodes(i) & x<=nodes(i+1);
    bin_num_x(in)=i;
end
sigma=interp_line(x,bin_num_x,GRID,param_nodes(:,1));
xi=interp_line(x,bin_num_x,GRID,param_nodes(:,2));

x=mod(x,360);
x(end)=360;

% plot quantiles
hold on; box on
plot(DATA.X,DATA.Y,'o','color',[1 1 1]*0.6)
plotcol(x,u,1,1,length(P)+2,'-',2)

Q=zeros(length(x),length(P));
for i=1:length(P)
    Q(:,i)=gpinv(P(i),xi,sigma,u);
    plotcol(x,Q(:,i),i+1,1,length(P)+2,'-',2)
end

yl=ylim;
for i=1:GRID.Nnode
    plot([1 1]*GRID.nodes.position(i),yl,'k--')
end
ylim(yl)
xlim([0 360])
xticks(0:90:360)

ylabel([DATA.name.Y ' [' DATA.unit.Y ']'])
xlabel([DATA.name.X ' [' DATA.unit.X ']'])
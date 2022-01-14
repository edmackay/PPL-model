function plot_grid_checks_1cov(DATA,GRID)

set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')

% shift gridded estimates for visualisation
ind=find(DATA.grid.x<GRID.nodes.position(1),1,'last');
if isempty(ind)
    ind=0;
end
x_grid=[DATA.grid.x(ind+1:end),DATA.grid.x(1:ind)+360];
sigma_grid=[DATA.grid.sigma(ind+1:end),DATA.grid.sigma(1:ind)];
xi_grid=[DATA.grid.xi(ind+1:end),DATA.grid.xi(1:ind)];

% check voronoi results against gridded estimates
figure
hold on; box on
plot(DATA.exceedance.X,DATA.exceedance.Z,'o','color',[1 1 1]*0.6)
plot(x_grid,sigma_grid,'r','linewidth',2)
plot(GRID.nodes.position,DATA.voronoi.xi_const.sigma([1:GRID.Nnode,1]),'k','linewidth',2)
plot(GRID.nodes.position,DATA.voronoi.xi_vary.sigma([1:GRID.Nnode,1]),'b','linewidth',2)
xlabel([DATA.name.X ' [' DATA.unit.X ']'])
ylabel('Threshold exceedance or GP scale')
xlim([GRID.nodes.position(1)-20 GRID.nodes.position(end)+20])
yl=ylim;
for i=1:GRID.Nnode+1
    h=plot([1 1]*GRID.nodes.position(i),yl,'k--');
    h.HandleVisibility='off';
end
ylim(yl)
L=legend('Threshold exceedance','$\sigma$ - local','Voronoi - constant shape','Voronoi - variable shape');
L.Box='off';

figure
hold on; box on
plot(x_grid,xi_grid,'r')
plot(GRID.nodes.position,DATA.voronoi.xi_const.xi+0*GRID.nodes.position,'k')
plot(GRID.nodes.position,DATA.voronoi.xi_vary.xi([1:GRID.Nnode,1]),'b')
xlabel([DATA.name.X ' [' DATA.unit.X ']'])
ylabel('GP shape')
xlim([GRID.nodes.position(1)-20 GRID.nodes.position(end)+20])
yl=ylim;
for i=1:GRID.Nnode+1
    h=plot([1 1]*GRID.nodes.position(i),yl,'k--');
    h.HandleVisibility='off';
end
ylim(yl)
legend('Local estimate','Voronoi - constant shape','Voronoi - variable shape')

end
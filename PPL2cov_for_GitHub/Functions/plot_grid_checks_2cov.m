function plot_grid_checks_2cov(DATA,GRID)

figure
hold on; grid on
scatter(DATA.exceedance.X(:,1),DATA.exceedance.X(:,2),'.')
triplot(GRID.tri.indices,GRID.nodes.position_repeated(:,1),GRID.nodes.position_repeated(:,2),'k')
scatter(GRID.nodes.position(:,1), GRID.nodes.position(:,2), 'ro', 'Filled')
xlabel([DATA.name.X1 ' [' DATA.unit.X1 ']'])
ylabel([DATA.name.X2 ' [' DATA.unit.X2 ']'])
axis square
xticks(-360:90:720)
yticks(-360:90:720)

% Triangulation projected onto margins
figure
subplot(1,2,1)
hold on; grid on; zoom on
scatter(DATA.exceedance.X(:,1),DATA.exceedance.Z)
yl=ylim;
for i=1:GRID.Nnode_unique
    plot([1 1]*GRID.nodes.position(i,1),yl,'k--')
end
xlabel([DATA.name.X1 ' [' DATA.unit.X1 ']'])
ylabel(['Threshold exceedance [' DATA.unit.Y ']'])

subplot(1,2,2)
hold on; grid on
scatter(DATA.exceedance.X(:,2),DATA.exceedance.Z)
yl=ylim;
for i=1:GRID.Nnode_unique
    plot([1 1]*GRID.nodes.position(i,2),yl,'k--')
end
xlabel([DATA.name.X2 ' [' DATA.unit.X2 ']'])
ylabel(['Threshold exceedance [' DATA.unit.Y ']'])

% 3D plot of GP scale estimates
figure
plot_initial_gridded_estimates_2cov_3D(DATA,'sigma',GRID)

% Overlay grid on initial GP scale estimates
figure
ax1=subplot(1,2,1);
plot_initial_gridded_estimates_2cov(DATA,'sigma')
triplot(GRID.tri.indices,GRID.nodes.position_repeated(:,1),GRID.nodes.position_repeated(:,2),'k')
title('GP scale - local')

ax2=subplot(1,2,2);
hold on
trisurf(GRID.tri.indices,GRID.nodes.position_repeated(:,1),GRID.nodes.position_repeated(:,2),0*GRID.nodes.number,DATA.voronoi.xi_const.sigma(GRID.nodes.number),'facecolor','interp')
view(2)
xlabel([DATA.name.X1 ' [' DATA.unit.X1 ']'])
ylabel([DATA.name.X2 ' [' DATA.unit.X2 ']'])
title('GP scale - Voronoi')
colorbar
axis([-360 720 -360 720])
axis square
zoom on
linkaxes([ax1 ax2])

figure
ax3=subplot(1,2,1);
plot_initial_gridded_estimates_2cov(DATA,'xi')
triplot(GRID.tri.indices,GRID.nodes.position_repeated(:,1),GRID.nodes.position_repeated(:,2),'k')
title('GP shape - local')
caxis([-0.5 0])

ax4=subplot(1,2,2);
hold on
trisurf(GRID.tri.indices,GRID.nodes.position_repeated(:,1),GRID.nodes.position_repeated(:,2),0*GRID.nodes.number,DATA.voronoi.xi_vary.xi(GRID.nodes.number),'facecolor','interp')
view(2)
xlabel([DATA.name.X1 ' [' DATA.unit.X1 ']'])
ylabel([DATA.name.X2 ' [' DATA.unit.X2 ']'])
title('GP shape - Voronoi')
colorbar
axis([-360 720 -360 720])
axis square
zoom on
linkaxes([ax3 ax4])
caxis([-0.5 0])
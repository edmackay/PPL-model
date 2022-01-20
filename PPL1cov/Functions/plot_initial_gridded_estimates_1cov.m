function plot_initial_gridded_estimates_1cov(DATA,type)

set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')

if strcmpi(type,'pdf')
    h=histogram(DATA.X,0:10:360,'Normalization','pdf');
    h.FaceColor=[1 1 1]*0.5;
    hold on; box on
    plot(DATA.grid.x,DATA.grid.pdf,'r','LineWidth',2)
    ylabel('Probability density');
    L=legend('Binned','Kernel density');
    L.Box='off';
elseif strcmpi(type,'thresh')
    plot(DATA.X,DATA.Y,'o','color',[1 1 1]*0.6)
    hold on; box on
    plot(DATA.grid.x,DATA.grid.thresh_raw,'r','linewidth',2)
    plot(DATA.grid.x,DATA.grid.thresh,'b','linewidth',2)
    ylabel([DATA.name.Y ' [' DATA.unit.Y ']'])
    L=legend('Data','Threshold - raw','Threshold - smoothed');
    L.Box='off';
elseif strcmpi(type,'sigma')
    plot(DATA.grid.x,DATA.grid.sigma)
    ylabel('GP scale')
elseif strcmpi(type,'xi')
    plot(DATA.grid.x,DATA.grid.xi)
    ylabel('GP shape')
end

xlabel([DATA.name.X ' [' DATA.unit.X ']'])
set(gca,'xtick',0:45:360)
xlim([0 360])

end
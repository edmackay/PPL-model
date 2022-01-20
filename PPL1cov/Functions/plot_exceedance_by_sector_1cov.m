function plot_exceedance_by_sector_1cov(DATA_obs,DATA_mod,type,xrange,yrange,FS)

set(0,'defaultTextInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex')

umax=max(DATA_obs.grid.thresh);
Xobs=DATA_obs.exceedance.X;
Yobs=DATA_obs.exceedance.Z+DATA_obs.exceedance.thresh;
Xmod=DATA_mod.X;
Ymod=DATA_mod.Y;
Xgrid=DATA_obs.grid.x;
Nboot=size(Ymod,2);

if nargin<6
    FS=10;
end

% sort data and calculate quantiles
nq=100;
Qi=logspace(-5,0,nq);
Ymod_Qi=zeros(nq,Nboot);
for i=1:Nboot
    if Nboot>1
        disp(['sorting bootstrap trial ' num2str(i)])
    end
    [Ymod(:,i),inds]=sort(Ymod(:,i));
    Xmod(:,i)=Xmod(inds,i);
    Ytemp=Ymod(Ymod(:,i)>umax,i);
    n=length(Ytemp);
    P=((1:n)-0.31)/(n+0.38);
    Q=(1-P);
    Ymod_Qi(:,i)=interp1(log(Q),Ytemp,log(Qi));
end

% crop to range where all simulations contain data
s=sum(isnan(Ymod_Qi),2)>1;
s0=find(s==0,1,'first');
s1=find(s==0,1,'last');
Qi=Qi(s0:s1);
Ymod_Qi=Ymod_Qi(s0:s1,:);

% calculate quantiles
Y025=quantile(Ymod_Qi,0.025,2);
Y50=quantile(Ymod_Qi,0.5,2);
Y975=quantile(Ymod_Qi,0.975,2);

% omnidirectional
if strcmpi(type(1:3),'sea')
    subplot(3,6,[1,2,7,8,13,14])
else
    subplot(3,5,[1,2,6,7,11,12])
end
hold on; box on; grid on
if Nboot==1
    exceedance_plot(Yobs(Yobs>umax),1,'ko')
    exceedance_plot(Ymod(Ymod>umax),0,'r-')
else
    exceedance_plot(Yobs(Yobs>umax),0,'ko')
    plot(Y50,Qi,'r-')
    h1=plot(Y025,Qi,'r--');
    h2=plot(Y975,Qi,'r--');
    h1.HandleVisibility='off';
    h2.HandleVisibility='off';
end
xlabel([DATA_obs.name.Y ' [' DATA_obs.unit.Y ']'])
title('All data')
L=legend('Observed','Model');
L.Box='off';
set(gca,'yMinorTick','off')
set(gca,'yminorgrid','off')
set(gca,'Fontsize',FS)
if nargin>4
    ylim(yrange)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(type(1:3),'dir')
    % directional
    Xobs(Xobs>337.5)=Xobs(Xobs>337.5)-360;
    Xmod(Xmod>337.5)=Xmod(Xmod>337.5)-360;
    Xgrid(Xgrid>337.5)=Xgrid(Xgrid>337.5)-360;
    
    b0=-22.5:45:315;
    sp=[4 5 10 15 14 13 8 3];
    str={'N','NE','E','SE','S','SW','W','NW'};
    for i=1:8
        bin_mod = Xmod>=b0(i) & Xmod<b0(i)+45;
        bin_obs = Xobs>=b0(i) & Xobs<b0(i)+45;
        bin_grid = Xgrid>=b0(i) & Xgrid<b0(i)+45;
        umax=max(DATA_obs.grid.thresh(bin_grid));
        
        Yobs_bin=Yobs(bin_obs);
        Ymod_bin=Ymod(bin_mod);
        
        subplot(3,5,sp(i))
        hold on; box on; grid on
        exceedance_plot(Yobs_bin(Yobs_bin>umax),0,'ko')
        exceedance_plot(Ymod_bin(Ymod_bin>umax),0,'r-')
        set(gca,'yMinorTick','off')
        set(gca,'yminorgrid','off')
        set(gca,'Fontsize',FS)
        ylabel('')
        title(str{i})
        if nargin>3
            xlim(xrange)
        end
        if nargin>4
            ylim(yrange)
        end
    end
    
    % plot the directional arrows in the centre
    subplot(3,5,9)
    hold on
    axis([-1 1 -1 1])
    set(gca,'dataaspectratio',[1 1 1],'plotboxaspectratiomode','auto')
    t=22.5:45:180;
    for i=1:4
        x0=cos(t(i)*pi/180);
        y0=sin(t(i)*pi/180);
        plot([-x0 x0],[-y0 y0],'k-')
    end
    r=0.95;
    text(0,r,'N','HorizontalAlignment','center','FontSize',FS)
    text(0,-r,'S','HorizontalAlignment','center','FontSize',FS)
    text(r,0,'E','HorizontalAlignment','center','FontSize',FS)
    text(-r,0,'W','HorizontalAlignment','center','FontSize',FS)
    r=0.65;
    text(r,r,'NE','HorizontalAlignment','center','FontSize',FS)
    text(r,-r,'SE','HorizontalAlignment','center','FontSize',FS)
    text(-r,r,'NW','HorizontalAlignment','center','FontSize',FS)
    text(-r,-r,'SW','HorizontalAlignment','center','FontSize',FS)
    axis off
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmpi(type(1:3),'sea')
    % seasonal
    b0=0:30:330;
    str={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
    sp=[3:6 9:12 15:18];
    for i=1:12
        bin_mod = Xmod>=b0(i) & Xmod<b0(i)+30;
        bin_obs = Xobs>=b0(i) & Xobs<b0(i)+30;
        bin_grid = Xgrid>=b0(i) & Xgrid<b0(i)+30;
        umax=max(DATA_obs.grid.thresh(bin_grid));
        
        Yobs_bin=Yobs(bin_obs);
        Ymod_bin=Ymod(bin_mod);
        
        subplot(3,6,sp(i))
        hold on; box on; grid on
        exceedance_plot(Yobs_bin(Yobs_bin>umax),0,'ko')
        exceedance_plot(Ymod_bin(Ymod_bin>umax),0,'r-')
        set(gca,'yMinorTick','off')
        set(gca,'yminorgrid','off')
        set(gca,'Fontsize',FS)
        ylabel('')
        title(str{i})
        if nargin>3
            xlim(xrange)
        end
        if nargin>4
            ylim(yrange)
        end
    end
end
function plot_exceedance_by_sector_2cov(DATA_obs,DATA_mod,xrange,yrange,FS,conf)

if nargin<5
    FS=10;
end
if nargin<6
    conf=false;
end

X1obs=DATA_obs.exceedance.X(:,1);
X2obs=DATA_obs.exceedance.X(:,2);
Yobs=DATA_obs.exceedance.Z+DATA_obs.exceedance.thresh;
X1mod=squeeze(DATA_mod.X(:,1,:));
X2mod=squeeze(DATA_mod.X(:,2,:));
Ymod=DATA_mod.Y;
Nboot=size(Ymod,2);

% sort data and calculate quantiles
nq=100;
Qi=logspace(-5,0,nq);
Ymod_Qi=zeros(nq,Nboot);
for i=1:Nboot
    if i>1
        disp(['Interpolating for trial ' num2str(i)])
    end
    [Ymod(:,i),inds]=sort(Ymod(:,i));
    X1mod(:,i)=X1mod(inds,i);
    X2mod(:,i)=X2mod(inds,i);
    n=length(Ymod(:,i));
    P=((1:n)-0.31)/(n+0.38);
    Q=(1-P);
    Ymod_Qi(:,i)=interp1(log(Q),Ymod(:,i),log(Qi));
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
figure
hold on; box on; grid on
if Nboot==1
    exceedance_plot(Yobs,0,'ko')
    exceedance_plot(Ymod,0,'r-')
else
    exceedance_plot(Yobs,0,'ko')
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
if nargin>3
    ylim(yrange)
end

% directional
X1obs(X1obs>337.5)=X1obs(X1obs>337.5)-360;
X1mod(X1mod>337.5)=X1mod(X1mod>337.5)-360;

b0=-22.5:45:315;
sp=[2 3 6 9 8 7 4 1];
str={'N','NE','E','SE','S','SW','W','NW'};
figure
if Nboot==1 || ~conf
    for i=1:8
        bin_mod = X1mod>=b0(i) & X1mod<b0(i)+45;
        bin_obs = X1obs>=b0(i) & X1obs<b0(i)+45;
        
        Yobs_bin=Yobs(bin_obs);
        Ymod_bin=Ymod(bin_mod);
        
        subplot(3,3,sp(i))
        hold on; box on; grid on
        exceedance_plot(Yobs_bin,0,'ko')
        exceedance_plot(Ymod_bin,0,'r-')
        set(gca,'yMinorTick','off')
        set(gca,'yminorgrid','off')
        set(gca,'Fontsize',FS)
        ylabel('')
        title(str{i})
        if nargin>2
            xlim(xrange)
        end
        if nargin>3
            ylim(yrange)
        end
    end
else
    for i=1:8
        bin_obs = X1obs>=b0(i) & X1obs<b0(i)+45;
        Yobs_bin=Yobs(bin_obs);
        
        nq=100;
        Qi=logspace(-5,0,nq);
        Ymod_Qi=zeros(nq,Nboot);
        for j=1:Nboot
            disp(['Interpolating for trial ' num2str(j)])
            bin_mod = X1mod(:,j)>=b0(i) & X1mod(:,j)<b0(i)+45;
            Ymod_bin=Ymod(bin_mod,j);
%             [Ymod(:,i),inds]=sort(Ymod(:,i));
%             X1mod(:,i)=X1mod(inds,i);
%             X2mod(:,i)=X2mod(inds,i);
            n=length(Ymod_bin);
            P=((1:n)-0.31)/(n+0.38);
            Q=(1-P);
            Ymod_Qi(:,j)=interp1(log(Q),Ymod_bin,log(Qi));
        end
        
        subplot(3,3,sp(i))
        hold on; box on; grid on
        exceedance_plot(Yobs_bin,0,'ko')
        plot(quantile(Ymod_Qi,0.5,2),Qi,'r-')
        plot(quantile(Ymod_Qi,0.025,2),Qi,'r--')
        plot(quantile(Ymod_Qi,0.975,2),Qi,'r--')
        set(gca,'yMinorTick','off')
        set(gca,'yminorgrid','off')
        set(gca,'Fontsize',FS)
        ylabel('')
        title(str{i})
        if nargin>2
            xlim(xrange)
        end
        if nargin>3
            ylim(yrange)
        end
    end
end

% plot the directional arrows in the centre
subplot(3,3,5)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% Seasonal exceedance

b0=0:30:330;
str={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
figure
for i=1:12
    bin_mod = X2mod>=b0(i) & X2mod<b0(i)+30;
    bin_obs = X2obs>=b0(i) & X2obs<b0(i)+30;
    
    Yobs_bin=Yobs(bin_obs);
    Ymod_bin=Ymod(bin_mod);
    
    subplot(3,4,i)
    hold on; box on; grid on
    exceedance_plot(Yobs_bin,0,'ko')
    exceedance_plot(Ymod_bin,0,'r-')
    set(gca,'yMinorTick','off')
    set(gca,'yminorgrid','off')
    set(gca,'Fontsize',FS)
    ylabel('')
    title(str{i})
    if nargin>2
        xlim(xrange)
    end
    if nargin>3
        ylim(yrange)
    end
end

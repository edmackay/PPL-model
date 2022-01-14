function [xout,Pout]=exceedance_plot(x,conf,ls,lw)

if nargin<2
    conf=0;
end
if nargin<3
    ls='ko';
end
if nargin<4
    lw=1;
end
% empirical exceedance probabilities
n=length(x);
x=sort(x);
k=1:n;
P=(k-0.31)/(n+0.38);

if nargout==0
    % plot results
    plot(x,1-P,ls,'LineWidth',lw)
    set(gca,'yscale','log')
    ylabel('Exceedance probability')
    % confidence bounds on exceedance probabilities
    if conf==1
        a=k;
        b=n-k+1;
        Plow=icdf('beta',0.025,a,b);
        Phigh=icdf('beta',0.975,a,b);
        
        hold on; box on; grid on
        h1=plot(x,1-Plow,'k--');
        h2=plot(x,1-Phigh,'k--');
        h1.HandleVisibility='off';
        h2.HandleVisibility='off';
%         legend('Data','95% conf','location','northeast')
    end
else
    xout=x;
    Pout=P;
end



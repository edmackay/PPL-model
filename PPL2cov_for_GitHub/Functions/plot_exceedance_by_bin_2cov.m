function plot_exceedance_by_bin_2cov(DATA_obs,DATA_mod,GRID,xrange,yrange)

% parse inputs
Ntri=GRID.Ntri_unique;
Uobs=DATA_obs.exceedance.thresh;
Umod=DATA_mod.U;
Yobs=DATA_obs.exceedance.Z+Uobs;
Ymod=DATA_mod.Y;
Nboot=size(Ymod,2);
bin_num_obs=DATA_obs.exceedance.bin_num;
bin_num_mod=DATA_mod.bin_num;

if Nboot==1
    for i=1:Ntri
        bin_obs = bin_num_obs==i;
        bin_mod = bin_num_mod==i;
        
        Yobs_bin=Yobs(bin_obs);
        Uobs_bin=Uobs(bin_obs);
        Ymod_bin=Ymod(bin_mod);
        Umod_bin=Umod(bin_mod);
        umax=max([Uobs_bin(:); Umod_bin(:)]);
        Yobs_bin=Yobs_bin(Yobs_bin>umax);
        Ymod_bin=Ymod_bin(Ymod_bin>umax);
        
        if numel(Yobs_bin)>5
            figure
            hold on; box on; grid on
            exceedance_plot(Yobs_bin,1,'ko')
            exceedance_plot(Ymod_bin,0,'r-')
            set(gca,'yminorgrid','off')
            set(gca,'yMinorTick','off')
            title(['Triangle ' num2str(i)])
            xlabel([DATA_obs.name.Y '[' DATA_obs.unit.Y ']'])
            if nargin>3
                xlim(xrange)
            end
            if nargin>4
                ylim(yrange)
            end
        end
    end
else
    for i=1:Ntri
        bin_obs = bin_num_obs==i;
        Yobs_bin=Yobs(bin_obs);
        Uobs_bin=Uobs(bin_obs);
        umax=max(Uobs_bin);
        Yobs_bin=Yobs_bin(Yobs_bin>umax);

        %%%%%%%%%%%%%%%
        % sort data and calculate quantiles
        nq=100;
        Qi=logspace(-5,0,nq);
        Ymod_Qi=zeros(nq,Nboot);
        for j=1:Nboot
            bin_mod = bin_num_mod(:,j)==i;
            Ymod_bin=Ymod(bin_mod,j);
            Ymod_bin=sort(Ymod_bin(Ymod_bin>umax));
            
            n=length(Ymod_bin);
            P=((1:n)-0.31)/(n+0.38);
            Q=(1-P);
            Ymod_Qi(:,j)=interp1(log(Q),Ymod_bin,log(Qi));
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
        %%%%%%%%%%%%%%%%%
        
        figure
        hold on; box on; grid on
        exceedance_plot(Yobs_bin,0,'ko')
        plot(Y50,Qi,'r')
        plot(Y025,Qi,'r--')
        plot(Y975,Qi,'r--')
        set(gca,'yminorgrid','off')
        set(gca,'yMinorTick','off')
        title(['Triangle ' num2str(i)])
        xlabel([DATA_obs.name.Y '[' DATA_obs.unit.Y ']'])
        if nargin>3
            xlim(xrange)
        end
        if nargin>4
            ylim(yrange)
        end
    end
end
function [sig, xi]=gp_fit_constrained(z,xi_min,xi_max)

% first guess is moment estimator of sigma and xi
zmax=max(z);
m=mean(z);
v=var(z);
sig0=m*(1-0.5*(1-m^2/v));
xi0=0.5*(1-m^2/v);
if xi0>xi_max
    xi0=xi_max;
end
if xi0<xi_min
    xi0=xi_min;
end
if xi0<0 && zmax>-sig0/xi0
    xi0=-sig0/zmax;
end
p0=[xi0, sig0];

% maximise likelihood starting from first guess
opts=optimset('display','off');
params=fminsearch(@(p)GPLikeConstrained(p,z,zmax,xi_min,xi_max),p0,opts);
% if exitflag<1
%     pause
% end
xi=params(1);
sig=params(2);

    % likelihood function
    function NLOGL=GPLikeConstrained(PARAMS,z,zmax,xi_min,xi_max)
        Xi=PARAMS(1);
        Sig=PARAMS(2);
        if Sig<0 || Xi<xi_min || Xi>xi_max || (Xi<0 && zmax>-Sig/Xi)
            NLOGL=Inf;
        else
            NLOGL = gplike([Xi,Sig],z);
        end
    end
end
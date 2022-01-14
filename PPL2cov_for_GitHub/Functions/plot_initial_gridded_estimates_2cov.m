function plot_initial_gridded_estimates_2cov(DATA,type,tight)

if nargin<3
    tight=false;
end

if strcmpi(type,'pdf')
    val=DATA.grid.pdf;
    str='Covariate PDF';
elseif strcmpi(type,'thresh')
    val=DATA.grid.thresh;
    str=['Threshold (smoothed) [' DATA.unit.Y ']'];
elseif strcmpi(type,'thresh_raw')
    val=DATA.grid.thresh_raw;
    str=['Threshold (unsmoothed) [' DATA.unit.Y ']'];    
elseif strcmpi(type,'sigma')
    val=DATA.grid.sigma;
    str='GP scale (smoothed)';
elseif strcmpi(type,'sigma_raw')
    val=DATA.grid.sigma_raw;
    str='GP scale (unsmoothed)';    
elseif strcmpi(type,'xi')
    val=DATA.grid.xi;
    str='GP shape';
end

% visualise gridded estimates
plot_covariate_image(DATA.grid.x1,DATA.grid.x2,val)
xlabel([DATA.name.X1 ' [' DATA.unit.X1 ']'])
ylabel([DATA.name.X2 ' [' DATA.unit.X2 ']'])
title(str)
colorbar
if tight
    axis([0 360 0 360])
    axis square
    pan on
else
    axis([-360 720 -360 720])
    axis square
    zoom on
end

    function plot_covariate_image(x,y,val)
        
        imagesc(x,y,val)
        hold on
        imagesc(x-360,y,val)
        imagesc(x+360,y,val)
        imagesc(x-360,y-360,val)
        imagesc(x+360,y-360,val)
        imagesc(x,y+360,val)
        imagesc(x,y-360,val)
        imagesc(x-360,y+360,val)
        imagesc(x+360,y+360,val)
        set(gca,'ydir','normal')
    end
end
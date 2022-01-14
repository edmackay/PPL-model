function DATA=set_initial_gridded_estimates_2cov(DATA)

% assumes periodic data in x1 and x2 with period 360

% initialise
x1=DATA.X(:,1);
x2=DATA.X(:,2);
ndat=length(x1);

% initialise output grid
nx1=length(DATA.grid.x1);
nx2=length(DATA.grid.x2);
THRESH_unsmooth =NaN(nx2,nx1);
PDF = NaN(nx2,nx1);

% calculate distance from x1 data to grid points
X1=repmat(x1,1,nx1);
X2=repmat(x2,1,nx1);
X1_grid=repmat(DATA.grid.x1,ndat,1);
DX1=abs(X1-X1_grid);
DX1(DX1>180)=360-DX1(DX1>180);
DX1square=DX1.^2;

bw=DATA.grid.KD_bw;
for i = 1:nx2
    disp(['Calculating rate and threshold for row ' num2str(i) '/' num2str(nx2)])

    % calculate distance from x2 data to grid points
    X2_grid=repmat(DATA.grid.x2(i),ndat,nx1);
    DX2=abs(X2-X2_grid);
    DX2(DX2>180)=360-DX2(DX2>180);
    distsquare=DX1square+DX2.^2;
    
    % find closest points and calculate quantile
    [~,inds]=sort(distsquare);
    inds=inds(1:DATA.grid.Nnearest,:);
    THRESH_unsmooth(i,:)=quantile(DATA.Y(inds),DATA.thresh_NEP);
    
    % kernel density estimate of PDF
    K=exp(-distsquare/(2*bw^2));
    PDF(i,:)=sum(K,1);
end
% normalise PDF
PDF=PDF/(ndat*2*pi*bw^2);

% smooth threshold using gaussian window
bw=DATA.grid.Thresh_binwidth;
THRESH_smooth=0*THRESH_unsmooth;
[X1_grid,X2_grid]=meshgrid(DATA.grid.x1,DATA.grid.x2);
for i = 1:nx2
    disp(['Smoothing threshold for row ' num2str(i) '/' num2str(nx2)])
    for j=1:nx1
        % calculate distance from x2 data to grid points
        DX1=abs(DATA.grid.x1(j)-X1_grid);
        DX1(DX1>180)=360-DX1(DX1>180);
        DX2=abs(DATA.grid.x2(i)-X2_grid);
        DX2(DX2>180)=360-DX2(DX2>180);
        distsquare=DX1.^2+DX2.^2;
        
        % smooth threshold using Gaussian window
        weight=exp(-distsquare/(2*bw^2));
        weight=weight/sum(weight(:));
        THRESH_smooth(i,j)=sum(sum(weight.*THRESH_unsmooth));
    end
end

% calculate threshold exceedances
U=interp2(DATA.grid.x1,DATA.grid.x2,THRESH_smooth,x1,x2);
inds=DATA.Y>U;

% set output structure
DATA.grid.thresh_raw=THRESH_unsmooth;
DATA.grid.thresh=THRESH_smooth;
DATA.grid.pdf=PDF;
DATA.exceedance.thresh=U(inds);
DATA.exceedance.Z=DATA.Y(inds)-DATA.exceedance.thresh;
DATA.exceedance.X=DATA.X(inds,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate GP parameters based on local exceedances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialise output grids
SIGMA_unsmooth =NaN(nx2,nx1);
XI =NaN(nx2,nx1);

% calculate distance from x1 data to grid points
x1=DATA.exceedance.X(:,1);
x2=DATA.exceedance.X(:,2);
ndat=length(x1);

X1=repmat(x1,1,nx1);
X2=repmat(x2,1,nx1);
X1_grid=repmat(DATA.grid.x1,ndat,1);
DX1=abs(X1-X1_grid);
DX1(DX1>180)=360-DX1(DX1>180);
DX1square=DX1.^2;

for i = 1:nx2
    disp(['Calculating GP parameters for row ' num2str(i) '/' num2str(nx2)])

    % calculate distance from x2 data to grid points
    X2_grid=repmat(DATA.grid.x2(i),ndat,nx1);
    DX2=abs(X2-X2_grid);
    DX2(DX2>180)=360-DX2(DX2>180);
    dist=DX1square+DX2.^2;
    
    % find closest points
    [~,inds]=sort(dist);
    inds=inds(1:DATA.grid.Nnearest,:);
    dataN=DATA.exceedance.Z(inds);
    
    % calculate moment estimators of GP params
    zmax=max(dataN);
    m=mean(dataN);
    v=var(dataN);
    xi=0.5*(1-m.^2./v);
    sigma=m.*(1-xi);
    bad=xi<0 & zmax>-sigma./xi;
    xi(bad)=-sigma(bad)./zmax(bad);

    SIGMA_unsmooth(i,:)=sigma;
    XI(i,:)=xi;
end

% smooth sigma using gaussian window
bw=DATA.grid.sigma_binwidth;
SIGMA_smooth=0*SIGMA_unsmooth;
[X1_grid,X2_grid]=meshgrid(DATA.grid.x1,DATA.grid.x2);
for i = 1:nx2
    disp(['Smoothing threshold for row ' num2str(i) '/' num2str(nx2)])
    for j=1:nx1
        % calculate distance from x2 data to grid points
        DX1=abs(DATA.grid.x1(j)-X1_grid);
        DX1(DX1>180)=360-DX1(DX1>180);
        DX2=abs(DATA.grid.x2(i)-X2_grid);
        DX2(DX2>180)=360-DX2(DX2>180);
        distsquare=DX1.^2+DX2.^2;
        
        % smooth threshold using Gaussian window
        weight=exp(-distsquare/(2*bw^2));
        weight=weight/sum(weight(:));
        SIGMA_smooth(i,j)=sum(sum(weight.*SIGMA_unsmooth));
    end
end

% update output structure
DATA.grid.sigma_raw=SIGMA_unsmooth;
DATA.grid.sigma=SIGMA_smooth;
DATA.grid.xi=XI;
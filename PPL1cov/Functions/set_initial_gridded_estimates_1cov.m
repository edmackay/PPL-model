function DATA=set_initial_gridded_estimates_1cov(DATA)

% assumes periodic data in x with period 360

% initialise
x=DATA.X;
ndat=length(x);

% initialise output grid
nx=length(DATA.grid.x);

% calculate distance from x1 data to grid points
X=repmat(x,1,nx);
X_grid=repmat(DATA.grid.x,ndat,1);
dist=abs(X-X_grid);
dist(dist>180)=360-dist(dist>180);

% kernel density estimate of PDF
bw=DATA.grid.KD_bw;
K=exp(-0.5*(dist/bw).^2);
PDF=sum(K,1)/(ndat*sqrt(2*pi)*bw);

% find closest points and calculate quantile
[~,inds]=sort(dist);
inds=inds(1:DATA.grid.Nnearest,:);
THRESH_unsmooth=quantile(DATA.Y(inds),DATA.thresh_NEP);

% smooth threshold using gaussian window
bw=DATA.grid.Thresh_bw;
THRESH_smooth=0*THRESH_unsmooth;
for i=1:length(DATA.grid.x)
	dist=abs(DATA.grid.x-DATA.grid.x(i));
    dist(dist>180)=360-dist(dist>180);
    weight=normpdf(dist,0,bw);
    weight=weight/sum(weight);
    THRESH_smooth(i)=sum(weight.*THRESH_unsmooth);
end

% calculate threshold exceedances
U=interp1(DATA.grid.x,THRESH_smooth,x);
inds=DATA.Y>U;

% set output structure
DATA.grid.thresh_raw=THRESH_unsmooth;
DATA.grid.thresh=THRESH_smooth;
DATA.grid.pdf=PDF;
DATA.exceedance.thresh=U(inds);
DATA.exceedance.Z=DATA.Y(inds)-U(inds);
DATA.exceedance.X=DATA.X(inds,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate GP parameters based on local exceedances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate distance from x1 data to grid points
x=DATA.exceedance.X(:,1);
ndat=length(x);

X=repmat(x,1,nx);
X_grid=repmat(DATA.grid.x,ndat,1);
dist=abs(X-X_grid);
dist(dist>180)=360-dist(dist>180);

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

% update output structure
DATA.grid.sigma=sigma;
DATA.grid.xi=xi;
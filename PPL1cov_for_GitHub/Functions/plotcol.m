function h=plotcol(x,y,col,cmin,cmax,lstyle,lwidth)

if numel(col)>1
    error('col must be a scalar')
end

if nargin<6
    lstyle='-';
end
if nargin<7
    lwidth=1;
end

cmap=colormap;
ncol=length(cmap);

% ncol=100;
% cmap(:,1)=linspace(1,0,ncol);
% cmap(:,2)=linspace(0,0,ncol);
% cmap(:,3)=linspace(0,1,ncol);
% colormap(cmap)

dc=(cmax-cmin)/(ncol-1);
line=round((col-cmin)/dc);
line(line<1)=1;
line(line>ncol)=ncol;
h1=plot(x,y,lstyle,'color',cmap(line,:),'linewidth',lwidth);
if nargout>0
    h=h1;
end

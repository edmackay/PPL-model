function plot_initial_gridded_estimates_2cov_3D(DATA,type,GRID,range)

if strcmpi(type,'pdf')
    val=DATA.grid.pdf;
    str='PDF';
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

% repeat gridded estimate
x=[DATA.grid.x1-360,DATA.grid.x1,DATA.grid.x1+360];
y=[DATA.grid.x2-360,DATA.grid.x2,DATA.grid.x2+360];
nx=length(DATA.grid.x1);
ny=length(DATA.grid.x2);
M=zeros(3*ny,3*nx);
for i=0:2
    for j=0:2
        M(i*ny+1:(i+1)*ny,j*nx+1:(j+1)*nx)=val;
    end
end

% find extends of GRID
if nargin<4
    xmin=min(GRID.nodes.position_repeated(:,1));
    xmax=max(GRID.nodes.position_repeated(:,1));
    ymin=min(GRID.nodes.position_repeated(:,2));
    ymax=max(GRID.nodes.position_repeated(:,2));
else
    xmin=range(1);
    xmax=range(2);
    ymin=range(3);
    ymax=range(4);
end

% crop grid
good=x>=xmin & x<=xmax;
x=x(good);
M=M(:,good);
good=y>=ymin & y<=ymax;
y=y(good);
M=M(good,:);

% interpolate gridded estimates onto triangulation
val_grid=interp2(DATA.grid.x1,DATA.grid.x2,val,mod(GRID.nodes.position(:,1),360),mod(GRID.nodes.position(:,2),360));

% visualise gridded estimates
surf(x,y,M,'EdgeColor','none');
hold on; rotate3d on
scatter3(GRID.nodes.position_repeated(:,1),GRID.nodes.position_repeated(:,2),val_grid(GRID.nodes.number),'ro','filled')
trisurf(GRID.tri.indices,GRID.nodes.position_repeated(:,1),GRID.nodes.position_repeated(:,2),val_grid(GRID.nodes.number),'facecolor','none')
xlabel([DATA.name.X1 ' [' DATA.unit.X1 ']'])
ylabel([DATA.name.X2 ' [' DATA.unit.X2 ']'])
zlabel(str)
light('Position',[0 0 3*max(DATA.Y)],'Style','local')
lighting gouraud
view(3)

end
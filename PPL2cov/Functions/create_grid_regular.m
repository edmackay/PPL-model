function [DATA, GRID]=create_grid_regular(DATA,marginal_nodes_x,marginal_nodes_y)

% Assumes periodic variables with period [0 360) for each variable
% nodes_in is an N x 2 array of x and y positions of nodes

% ensure nodes and data are on [0 360)
marginal_nodes_x=sort(mod(marginal_nodes_x,360));
marginal_nodes_y=sort(mod(marginal_nodes_y,360));

% shift data to start at first node on each axis
X=mod(DATA.exceedance.X,360);
shiftx=X(:,1)<marginal_nodes_x(1);
shifty=X(:,2)<marginal_nodes_y(1);
X(shiftx,1)=X(shiftx,1)+360;
X(shifty,2)=X(shifty,2)+360;

% ensure user-defined nodes are sorted column vectors
marginal_nodes_x=sort(marginal_nodes_x(:));
marginal_nodes_y=sort(marginal_nodes_y(:));

% count number of nodes in each covariate
nx=length(marginal_nodes_x);
ny=length(marginal_nodes_y);

% repeat first values
marginal_nodes_x=[marginal_nodes_x; marginal_nodes_x(1)+360];
marginal_nodes_y=[marginal_nodes_y; marginal_nodes_y(1)+360];

% find mid points
xmid=0.5*(marginal_nodes_x(1:end-1)+marginal_nodes_x(2:end));
ymid=0.5*(marginal_nodes_y(1:end-1)+marginal_nodes_y(2:end));

% create grid of node positions
[nodex_grid, nodey_grid]=meshgrid(marginal_nodes_x,marginal_nodes_y);
[nodexmid_grid, nodeymid_grid]=meshgrid(xmid,ymid);
nodes_out=[nodex_grid(:), nodey_grid(:);
    nodexmid_grid(:), nodeymid_grid(:)];

% create list of triangle vertices
ntri=nx*ny*4;
TRI=zeros(ntri,3);
ind=0;
for i=1:nx
    for j=1:ny
        ind1=find(nodes_out(:,1)==marginal_nodes_x(i) & nodes_out(:,2)==marginal_nodes_y(j));
        ind2=find(nodes_out(:,1)==marginal_nodes_x(i+1) & nodes_out(:,2)==marginal_nodes_y(j));
        ind3=find(nodes_out(:,1)==marginal_nodes_x(i+1) & nodes_out(:,2)==marginal_nodes_y(j+1));
        ind4=find(nodes_out(:,1)==marginal_nodes_x(i) & nodes_out(:,2)==marginal_nodes_y(j+1));
        indmid=find(nodes_out(:,1)==xmid(i) & nodes_out(:,2)==ymid(j));
        
        TRI(ind+1,:)=[ind1 ind2 indmid];
        TRI(ind+2,:)=[ind2 ind3 indmid];
        TRI(ind+3,:)=[ind3 ind4 indmid];
        TRI(ind+4,:)=[ind4 ind1 indmid];
        ind=ind+4;
    end
end

% set node numbers
repeated = nodes_out(:,1)==max(nodes_out(:,1)) | nodes_out(:,2)==max(nodes_out(:,2));
in = ~repeated;
nodes_in=nodes_out(in,:);
node_num=zeros(length(nodes_out),1);
node_num(in)=1:sum(in);
repeat=find(~in);
for i=1:length(repeat)
    pos=mod(nodes_out(repeat(i),:),360);
    node_num(repeat(i))=find(nodes_in(:,1)==pos(1) & nodes_in(:,2)==pos(2));
end

% figure
% hold on
% scatter(X(:,1),X(:,2),'.')
% triplot(TRI,nodes_out(:,1),nodes_out(:,2),'k')
% scatter(nodes_in(:,1),nodes_in(:,2),'ko','filled')
% for i=1:length(nodes_out)
%     text(nodes_out(i,1)+5,nodes_out(i,2),num2str(node_num(i)),'color','r','FontSize',14)
% end

% calculate interpolation coefficients and bin numbers
bin_num=zeros(length(X),1);
ntri=length(TRI);
interp_coef=zeros(3,3,ntri);
for i=1:ntri
    M=[ones(3,1) nodes_out(TRI(i,:),:)];
    interp_coef(:,:,i)=inv(M);
    in = inpolygon(X(:,1),X(:,2),nodes_out(TRI(i,:),1),nodes_out(TRI(i,:),2));
    bin_num(in)=i;
end

% calculate distance to nodes
dist=zeros(length(X),length(nodes_in));
for i=1:length(nodes_in)
    dx=abs(X(:,1)-nodes_in(i,1));
    dy=abs(X(:,2)-nodes_in(i,2));
    dx(dx>180)=360-dx(dx>180);
    dy(dy>180)=360-dy(dy>180);
    dist(:,i)=dx.^2+dy.^2;
end
[~,nearest_node]=min(dist,[],2);

% update data structure
DATA.exceedance.X=X;
DATA.exceedance.bin_num=bin_num;
DATA.exceedance.bin_num_unique=bin_num;
DATA.exceedance.nearest_node=nearest_node;
DATA.exceedance.Ndat=length(X);

% create grid structure
GRID.isregular=true;
GRID.tri.indices=TRI;
GRID.tri.number=(1:ntri)';
GRID.nodes.position=nodes_in;
GRID.nodes.position_repeated=nodes_out;
GRID.nodes.number=node_num;
GRID.interp_coef=interp_coef;
GRID.Ntri_unique=ntri;
GRID.Ntri_repeated=ntri;
GRID.Nnode_unique=length(nodes_in);
GRID.Nnode_repeated=length(nodes_out);

end
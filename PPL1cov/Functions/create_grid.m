function [DATA, GRID]=create_grid(DATA,nodes)

% Assumes periodic covariate with period [0 360)

% ensure nodes and data are on [0 360)
nodes=sort(mod(nodes(:),360));

% shift data to start at first node
X=mod(DATA.exceedance.X,360);
shift=X<nodes(1);
X(shift)=X(shift)+360;

% repeat first value
Nnode=length(nodes);
nodes=[nodes; nodes(1)+360];
node_num=[(1:Nnode)'; 1];

% calculate interpolation coefficients and bin numbers
bin_num=zeros(length(X),1);
for i=1:Nnode
    in = X>=nodes(i) & X<=nodes(i+1);
    bin_num(in)=i;
end

% calculate distance to nodes
dist=zeros(length(X),Nnode);
for i=1:Nnode
    dx=abs(X-nodes(i));
    dx(dx>180)=360-dx(dx>180);
    dist(:,i)=dx.^2;
end
[~,nearest_node]=min(dist,[],2);

% update data structure
DATA.exceedance.X=X;
DATA.exceedance.bin_num=bin_num;
DATA.exceedance.nearest_node=nearest_node;
DATA.exceedance.Ndat=length(X);

% create grid structure
GRID.nodes.position=nodes;
GRID.nodes.number=node_num;
GRID.Nnode=Nnode;

end
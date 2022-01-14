function [DATA, GRID]=create_grid_delaunay(DATA,nodes_in)

% Assumes periodic variables with period [0 360) for each variable
% nodes_in is an N x 2 array of x and y positions of nodes

% ensure nodes and data are on [0 360)
nodes_in=mod(nodes_in,360);
X=mod(DATA.exceedance.X,360);

% repeat nodes shifted by +/- 360
N=length(nodes_in);
node_num_initial=(1:N)';
node_num_initial=repmat(node_num_initial,9,1);
nodes_repeated=repmat(nodes_in,9,1);
px=[-1 0 1];
py=[-1 0 1];
[px,py]=meshgrid(px,py);
px=px(:);
py=py(:);
for i=1:9
    nodes_repeated((i-1)*N+1:i*N,1)=nodes_repeated((i-1)*N+1:i*N,1)+px(i)*360;
    nodes_repeated((i-1)*N+1:i*N,2)=nodes_repeated((i-1)*N+1:i*N,2)+py(i)*360;
end

% create triangulation
TRI_initial=delaunay(nodes_repeated(:,1),nodes_repeated(:,2));

% find out which triangles contain data
nbin=size(TRI_initial,1);
notempty=false(nbin,1);
for i=1:nbin
    in = inpolygon(X(:,1),X(:,2),nodes_repeated(TRI_initial(i,:),1),nodes_repeated(TRI_initial(i,:),2));
    if sum(in)>0
        notempty(i)=1;
    end
end
TRI_notempty=TRI_initial(notempty,:);

% crop to only nodes that are used
inds=unique(TRI_notempty(:));
nodes_out=nodes_repeated(inds,:);
node_num=node_num_initial(inds);
TRI=TRI_notempty;
for i=1:length(inds)
    j=TRI_notempty==inds(i);
    TRI(j)=i;
end
ntri=length(TRI);

% reorder triangles to ensure numbering is counter-clockwise from bottom left
for i=1:length(TRI)
    % get vertex positions
    vertx=nodes_out(TRI(i,:),1)';
    verty=nodes_out(TRI(i,:),2)';
    % find index of lower left corner
    rho=vertx+verty;
    ind=find(rho==min(rho));
    if length(ind)==2
       [~,indx]=min(vertx(ind));
       ind=ind(indx);
    end
    % reorder if necessary
    if ind==2
        TRI(i,:)=TRI(i,[2 3 1]);
        vertx=vertx([2 3 1]);
        verty=verty([2 3 1]);
    elseif ind==3
        TRI(i,:)=TRI(i,[3 1 2]);
        vertx=vertx([3 1 2]);
        verty=verty([3 1 2]);
    end
    % check if triangle vertices are order counter-clockwise
    AB=[vertx(2)-vertx(1), verty(2)-verty(1), 0];
    AC=[vertx(3)-vertx(1), verty(3)-verty(1), 0];
    c=cross(AB,AC);
    if c(3)<0
        TRI(i,:)=TRI(i,[1 3 2]);
        vertx=vertx([1 3 2]);
        verty=verty([1 3 2]);
    end
    
%     figure
%     hold on; grid on
%     scatter(vertx,verty)
%     text(vertx(1)+5,verty(1),'1')
%     text(vertx(2)+5,verty(2),'2')
%     text(vertx(3)+5,verty(3),'3')
%     axis([-360 720 -360 720])
%     axis_equal_and_square
end

% find unique triangles
TRI_number=NaN(ntri,1);
unique_vertices=NaN(ntri,6);
unique_vertices(1,:)=[nodes_out(TRI(1,:),1)' nodes_out(TRI(1,:),2)'];
unique_tri_count=1;
TRI_number(1)=unique_tri_count;
for i=2:ntri
    % create shifted copies of triangle
    shift=nan(9,6);
    for j=1:9
        shift(j,:)=[nodes_out(TRI(i,:),1)'+px(j)*360, nodes_out(TRI(i,:),2)'+py(j)*360];
    end
    % see if triangle is shifted copy of previous triangles
    [a,b]=ismember(shift,unique_vertices(1:unique_tri_count,:),'rows');
    
%     figure
%     hold on
%     triplot(TRI(i,:),nodes_out(:,1),nodes_out(:,2),'k')
%     for j=1:9
%         plot(shift(j,[1 2 3 1]),shift(j,[4 5 6 4]),'b')
%     end
%     for j=1:count
%         plot(unique_vertices(j,[1 2 3 1]),unique_vertices(j,[4 5 6 4]),'r--')
%         cx=sum(unique_vertices(j,1:3))/3;
%         cy=sum(unique_vertices(j,4:6))/3;
%         text(cx,cy,num2str(j),'Color','r','FontSize',14)
%     end
    
    if sum(a)==0
        unique_tri_count=unique_tri_count+1;
        TRI_number(i)=unique_tri_count;
        unique_vertices(unique_tri_count,:)=[nodes_out(TRI(i,:),1)' nodes_out(TRI(i,:),2)'];
    else
        TRI_number(i)=b(find(b));
    end
end

% figure
% hold on
% % grid on
% set(gca,'xtick',-360:180:720)
% set(gca,'ytick',-360:180:720)
% axis([-360 720 -360 720])
% axis square
% % plot([0 360 360 0 0],[0 0 360 360 0],'k')
% scatter(X(:,1),X(:,2),'.')
% scatter(nodes_in(:,1),nodes_in(:,2),'ko','filled')
% triplot(TRI_initial,nodes_repeated(:,1),nodes_repeated(:,2),'k')
% h=trisurf(TRI_notempty,nodes_repeated(:,1),nodes_repeated(:,2),0*nodes_repeated(:,2));
% h.EdgeColor='none';
% h.FaceAlpha=0.5;
% for i=1:length(TRI)
%     cx=sum(nodes_out(TRI(i,:),1))/3;
%     cy=sum(nodes_out(TRI(i,:),2))/3;
%     text(cx,cy,num2str(TRI_number(i)),'Color','r','FontSize',14)
% end

% calculate interpolation coefficients and bin numbers
bin_num=zeros(length(X),1);
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
DATA.exceedance.bin_num_unique=TRI_number(bin_num);
DATA.exceedance.nearest_node=nearest_node;
DATA.exceedance.Ndat=length(X);

% create grid structure
GRID.isregular=false;
GRID.tri.indices=TRI;
GRID.tri.number=TRI_number;
GRID.nodes.position=nodes_in;
GRID.nodes.position_repeated=nodes_out;
GRID.nodes.number=node_num;
GRID.interp_coef=interp_coef;
GRID.Ntri_unique=max(TRI_number);
GRID.Ntri_repeated=length(TRI);
GRID.Nnode_unique=length(nodes_in);
GRID.Nnode_repeated=length(nodes_out);

end
function [val_X,a]=interp_tri(X,bin_num,GRID,val_node)

% set values at repeated nodes
val_node_all=val_node(GRID.nodes.number);

% calculate for each triangle
val_X = NaN(length(X),1);
a = NaN(GRID.Ntri_repeated,3);
for i = 1:GRID.Ntri_repeated
    % calculate weights
    a(i,:) = GRID.interp_coef(:,:,i)*val_node_all(GRID.tri.indices(i,:));
    % interpolate within bin
    inds = bin_num==i;
    val_X(inds) = a(i,1) + X(inds,:)*(a(i,2:end))';
end


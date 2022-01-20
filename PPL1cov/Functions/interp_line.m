function [val_X,a]=interp_line(X,bin_num,GRID,val_node)

% relative distance from first node
d=(X-GRID.nodes.position(bin_num))./(GRID.nodes.position(bin_num+1)-GRID.nodes.position(bin_num));

% repeat first node value (assuming period covariate)
val_node=[val_node(:); val_node(1)];

% interpolating value
val_X=(1-d).*val_node(bin_num) + d.*val_node(bin_num+1);

% gradients in each bin
a=diff(val_node(GRID.nodes.number))./diff(GRID.nodes.position);
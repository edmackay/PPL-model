function [DATA, GRID, filename]=create_grid(DATA,reg_grid,nodes)

if reg_grid
    [DATA, GRID]=create_grid_regular(DATA,nodes.x1,nodes.x2);
    filename=[DATA.name.dataset '_binned_regular_' num2str(length(nodes.x1)) 'x' num2str(length(nodes.x2))];
else
    [DATA, GRID]=create_grid_delaunay(DATA,nodes);
    filename=[DATA.name.dataset '_binned_irregular_' num2str(size(nodes,1))];
end



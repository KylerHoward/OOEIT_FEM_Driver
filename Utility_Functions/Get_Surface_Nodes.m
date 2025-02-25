function [surf_nodes, surf_ind] = Get_Surface_Nodes(G_nodes, surf_connects)
    %{
    Extract just the surface nodes from a given surface mesh
    10/4/24 - Kyler Howard

    param: G_nodes       - All nodes in the global mesh
    param: surf_connects - Connectivity matrix for the surface faces

    return: surf_nodes - Matrix of the x/y/z coordinates for surface nodes
    return: surf_ind   - Indicies of the global nodes that are on the surface
    %}

    % Initalize nodes for the trachea
    surf_nodes = zeros(size(unique(surf_connects),1), 3);

    i = 1;
    % Loop through the connectivity matrix to find node coords
    surf_ind = unique(surf_connects, "stable");
    for node = surf_ind'
        surf_nodes(i,:) = G_nodes(node,:);
        i                = i + 1;
    end
end

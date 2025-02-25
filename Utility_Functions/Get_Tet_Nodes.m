function [tissue_nodes, node_ind] = Get_Tet_Nodes(G_nodes, tissue_connect)
    %{
    Go through and find all faces on the original mesh that contain the electrode nodes
    10/4/24 - Kyler Howard

    param: G_nodes        - All nodes in the global mesh
    param: tissue_connect - Connectivity matrix for an individual tissue

    return: tissue_nodes - Matrix of the x/y/z coordinates for tissue nodes
    return: node_ind     - Indicies of the global nodes that are in the tissue
    %}

    % Initalize nodes for the tissue
    tissue_nodes = zeros(size(unique(tissue_connect),1), 3);
    
    i = 1;
    % Loop through the connectivity matrix to find node coords
    node_ind = unique(tissue_connect, "stable");
    for node = node_ind'
        tissue_nodes(i,:) = G_nodes(node,:);
        i                 = i + 1;
    end
end
function E_connect = Align_Electrode_Faces2(G_nodes, surface_faces, E_nodes)
    %{
    Go through and find all faces on the original mesh that contain the electrode nodes
    10/4/24 - Kyler Howard

    param: G_nodes - All nodes in the global mesh
    param: surface_faces - Just the surface faces of the soft tissue label
    param: E_nodes - Array containing the nodes for A electrode
    param: flags - Various flags controlling plotting and other parameters

    return: E_connect - Cell array containing the connectivity matrix for each electrode
    %}
    
    % Find the global nodes indicies for the electrode
    [~, indices] = intersect(G_nodes, E_nodes, 'rows', 'stable');
    
    % Determine if the surface faces contain any of the indices
    isInIndices = ismember(surface_faces, indices);
    
    % Sum the logical matrix along the columns (how many of the indicies are in each row)
    rowSum = sum(isInIndices, 2);
    
    % Extract the rows which have all three indicies in them
    filteredRows = surface_faces(rowSum == 3, :);
    E_connect    = filteredRows;
end
function E_connect = Align_Electrode_Faces(G_nodes, surface_faces, E_nodes, flags)
    %{
    Go through and find all faces on the original mesh that contain the electrode nodes
    10/4/24 - Kyler Howard

    param: G_nodes - All nodes in the global mesh
    param: surface_faces - Just the surface faces of the soft tissue label
    param: E_nodes - Cell array containing the nodes for each electrode
    param: flags - Various flags controlling plotting and other parameters

    return: E_connect - Cell array containing the connectivity matrix for each electrode
    %}

    % Initialize the connectivity cell array
    E_connect = cell(size(E_nodes));
    
    % Prepare a figure if plotting electrodes
    if flags.plot_electrodes == 1
        figure
        hold on
        trimesh(surface_faces, G_nodes(:,1),G_nodes(:,2),G_nodes(:,3))
    end
    
    % Loop through each electrode
    for l = 1:length(E_nodes)
    
        % Find the global nodes indicies for the electrode
        [~, indices] = intersect(G_nodes, E_nodes{l}, 'rows', 'stable');
        
        % Determine if the surface faces contain any of the indices
        isInIndices = ismember(surface_faces, indices);
        
        % Sum the logical matrix along the columns (how many of the indicies are in each row)
        rowSum = sum(isInIndices, 2);
        
        % Extract the rows which have all three indicies in them
        filteredRows = surface_faces(rowSum == 3, :);
        E_connect{l} = filteredRows;
        
        % Plot the individual electrode if plotting electrodes
        if flags.plot_electrodes == 1
            scatter3(E_nodes{l}(:,1), E_nodes{l}(:,2), E_nodes{l}(:,3),'r', 'filled')
            trimesh(filteredRows, G_nodes(:,1),G_nodes(:,2),G_nodes(:,3), 'FaceColor', [0.6,0,0], 'EdgeColor', [0.6,0,0])
            xlabel("X (mm)");
            ylabel("Y (mm)");
            zlabel("Z (mm)");
            axis equal
        end
    end
end
function E_connect = Align_Electrode_Faces(G_nodes, surface_faces, E_nodes, flags)
    %{
    Go through and find all faces on the original mesh that contain the electrode nodes
    Plot the electrodes with their face normals
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
        trimesh(surface_faces, G_nodes(:,1),G_nodes(:,2),G_nodes(:,3),"FaceColor","cyan","EdgeColor",'blue')
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
            % Find center and normals for each electrode & plot them
            [centers, normals] = Find_Plot_Normals(E_connect{l}, G_nodes);
            
            % Check if the normals are pointing inward (>0), or outward (<0)
            if dot(mean(normals), mean(G_nodes) - mean(centers)) < 0
                scatter_color = [1,  0,0];
                mesh_color    = [0.6,0,0];
            else
                scatter_color = [0,1,  0];
                mesh_color    = [0,0.6,0];
            end
            
            scatter3(E_nodes{l}(:,1), E_nodes{l}(:,2), E_nodes{l}(:,3), 'MarkerEdgeColor', scatter_color, 'MarkerFaceColor', scatter_color)
            trimesh(filteredRows, G_nodes(:,1),G_nodes(:,2),G_nodes(:,3), 'FaceColor', mesh_color, 'EdgeColor', mesh_color)
            % quiver3(centers(:,1),centers(:,2),centers(:,3), normals(:,1),normals(:,2),normals(:,3), 2, 'color',scatter_color);
            xlabel("X (mm)");
            ylabel("Y (mm)");
            zlabel("Z (mm)");
            axis equal
            title("Red faces are outward normals")
        end
    end
end
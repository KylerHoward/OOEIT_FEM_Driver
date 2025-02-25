function [E_nodes, E_connect, perim_mm] = Make_Electrodes(boundary_nodes, carina_height, flags)
    %{
    Find the nodes that make up all 32 electrodes, for either a belt or patch configuration
    9/27/24 - Kyler Howard

    param: boundary_nodes - All nodes on the boundary of the surface
    param: carina_height  - Height of the carina for reference of electrode placement
    param: flags          - Various flags controlling plotting and other parameters

    return: E_nodes   - 1x32 Cell array containing electrode nodes
    return: E_connect - 1x32 Cell array containing Delauny electrode triangulation
    %}

% clear
% clc
% close all
% load electrode_temp.mat 
% flags.E_choice = 3;
% ----------------------------------------------------------------------- %
%% ------------------------------- Setup -------------------------------- %
% ----------------------------------------------------------------------- %
    % 1: Large patch front back
    % 2: Small patch front back
    % 3: Two rows of large belts
    % 4: Two rows of small belts

    % Constructing each electrode type
    L_square.type       = "patch";
    L_square.E_width    = 10;               % mm
    L_square.E_height   = 10;               % mm
    L_square.E_count    = [4,4];            % Electrodes per row and per column
    L_square.gap_width  = 2.5;              % mm (edge to edge)
    L_square.gap_height = 2.5;              % mm (edge to edge)
    
    S_square.type       = "patch";
    S_square.E_width    = 7;                % mm
    S_square.E_height   = 7;                % mm
    S_square.E_count    = [4,4];            % Electrodes per row and per column
    S_square.gap_width  = 2.5;              % mm (edge to edge)
    S_square.gap_height = 2.5;              % mm (edge to edge)
    
    L_belt.type    = "belt";
    L_belt.E_dia   = 15;                    % mm
    L_belt.E_rad   = L_belt.E_dia / 2;      % mm
    L_belt.E_count = 8;                     % Electrodes per belt
    L_belt.gap     = 5;                     % mm (edge to edge) 
    
    S_belt.type    = "belt";
    S_belt.E_dia   = 12;                    % mm
    S_belt.E_rad   = S_belt.E_dia / 2;      % mm
    S_belt.E_count = 8;                     % Electrodes per belt
    S_belt.gap     = 5;                     % mm (edge to edge)

    choices = {L_square, S_square, L_belt, S_belt};
    E = choices{flags.E_choice};
    
% ----------------------------------------------------------------------- %
%% --------------------------- Center Nodes ----------------------------- %
% ----------------------------------------------------------------------- %
    % Find the plane in which the center lies
    try 
        plane_gap = E.E_height / 2;
    catch
        plane_gap = E.E_rad;
    end
    z_check = (boundary_nodes(:,3) > carina_height - plane_gap) & (boundary_nodes(:,3) < carina_height + plane_gap);
    plane   = boundary_nodes(z_check,:);

    if E.type == "belt"
        height_low  = carina_height - E.E_dia - E.gap;
        z_check_low = (boundary_nodes(:,3) > height_low - plane_gap) & (boundary_nodes(:,3) < height_low + plane_gap);
        plane_low   = boundary_nodes(z_check_low,:);
    end

    % Find the nodes that are the center of the front/back of the body
    x_mid = (min(plane(:,1)) + max(plane(:,1))) / 2;
    y_min = min(plane(:,2));
    y_max = max(plane(:,2));
    front = plane(dsearchn(plane, [x_mid, y_min, carina_height]),:);
    back  = plane(dsearchn(plane, [x_mid, y_max, carina_height]),:);
    % front = find_node(plane, x_mid, y_min, carina_height);
    % back  = find_node(plane, x_mid, y_max, carina_height);

% ----------------------------------------------------------------------- %
%% ------------------------- Finding Perimeter -------------------------- %
% ----------------------------------------------------------------------- %
    center = (max(plane,[],1) + min(plane,[],1)) / 2;

    i = 1;
    point = zeros(100, 3);
    perim_mm = 0;
    for theta = 0 : (2*pi)/100 : 2*pi - (2*pi)/100 
        radius =  Parameratize_Bdry(plane, 15, theta);
        point(i,:) = [center(1) + radius*cos(theta), center(2) + radius*sin(theta), center(3)];
        
        if i >= 2  
            dist = sqrt((point(i-1,1) - point(i,1))^2 + (point(i-1,2) - point(i,2))^2);
            perim_mm = perim_mm + dist;
        end

        i = i + 1;
    end

% ----------------------------------------------------------------------- %
%% ---------------------- Finding Electrode Nodes ----------------------- %
% ----------------------------------------------------------------------- %
    if flags.E_choice <= 2 % PATCHES
        % Find the corners of the electrodes
        [front_nodes, front_connect] = create_patch(boundary_nodes, front, E, "front");
        [back_nodes,  back_connect]  = create_patch(boundary_nodes, back,  E, "back");
        E_nodes                      = [front_nodes; back_nodes];
        E_connect                    = [front_connect; back_connect];

    else % BELTS
        [top_nodes, top_connect] = create_belt(boundary_nodes, plane, E);
        [bot_nodes, bot_connect] = create_belt(boundary_nodes, plane_low, E);
        E_nodes                  = cat(1, bot_nodes,   top_nodes)';
        E_connect                = cat(1, bot_connect, top_connect)';
    end

% ----------------------------------------------------------------------- %
%% ------------------------------ Plotting ------------------------------ %
% ----------------------------------------------------------------------- %

    if flags.plot_electrodes == 1 && E.type == "patch"
        figure(); 
            hold on; 
            for cell_i = 1:size(front_nodes,1)
                scatter3(front_nodes{cell_i}(:,1), front_nodes{cell_i}(:,2), front_nodes{cell_i}(:,3), 'filled');
            end
            % scatter3(back_nodes(:,1), back_nodes(:,2), back_nodes(:,3), 'r', 'filled')
            for cell_i = 1:size(back_nodes,1)
                scatter3(back_nodes{cell_i}(:,1), back_nodes{cell_i}(:,2), back_nodes{cell_i}(:,3), 'filled');
            end
            legend('Location','eastoutside')
            xlabel('X (mm)')
            ylabel('Y (mm)')
            zlabel('Z (mm)')
            axis equal
    elseif flags.plot_electrodes == 1 && E.type == "belt"
        figure(); 
            hold on
            for node_i = 1:length(E_nodes)
                scatter3(E_nodes{node_i}(:,1), E_nodes{node_i}(:,2), E_nodes{node_i}(:,3), 'filled')
            end
            legend('Location','eastoutside')
            xlabel('X (mm)')
            ylabel('Y (mm)')
            zlabel('Z (mm)')
            axis equal
    end
end
% ----------------------------------------------------------------------- %
%% -------------------------- Custom Functions -------------------------- %
% ----------------------------------------------------------------------- %

function [center, plane] = find_center(nodes, height, FoB)
    %{
    Find the center of the body and the plane it is on
    9/27/24 - Kyler Howard

    param: nodes  - All nodes on the boundary of the surface
    param: height - The wanted height for the plane
    param: FoB    - Flag for if we are looking at the front or back of the body

    return: center - Coordinate of the center
    return: plane  - Plane containing the center
    %}

    z_check = (nodes(:,3) > height - 1.5) & (nodes(:,3) < height + 1.5);
    plane   = nodes(z_check,:);

    % Find the nodes that are the center of the front/back of the body
    x_mid = (min(plane(:,1)) + max(plane(:,1))) / 2;
    y_min = min(plane(:,2));
    y_max = max(plane(:,2));
    if lower(FoB) == "front"
        center = plane(dsearchn(plane, [x_mid, y_min, height]),:);
    elseif lower(FoB) == "back"
        center  = plane(dsearchn(plane, [x_mid, y_max, height]),:);
    end
end

function [E_nodes, E_connect] = create_electrode(local_nodes, coord, E)
    %{
    Create an individual electrode and find the nodes/connectivity
    9/27/24 - Kyler Howard

    param: local_nodes - All nodes on to look at for making the electrode
    param: coord       - Reference coordinate for each electrode
                         Bottom left corner  for patches
                         Center of electrode for belt
    param: E           - Structure containing the size and type of electrode

    return: corners   - Coordinates of the electrode
    return: E_connect - Delauny trainagulation for the electrode
    %}

    x = coord(1); y = coord(2); z = coord(3);
    if E.type == "patch" % PATCH ELECTRODE
        % Find the corners of one electrode
        flat_corners = [x + E.gap_width/2,             y, z + E.E_height/2;...
                        x + E.gap_width/2 + E.E_width, y, z + E.E_height/2;...
                        x + E.gap_width/2,             y, z - E.E_height/2;...
                        x + E.gap_width/2 + E.E_width, y, z - E.E_height/2];

        E_nodes      = local_nodes(dsearchn(local_nodes, flat_corners),:);

        xs  = local_nodes(:,1); 
        ys  = local_nodes(:,2); 
        zs  = local_nodes(:,3);
        ind = xs >= min(E_nodes(:,1)) & xs <= max(E_nodes(:,1)) & ys >= min(E_nodes(:,2)) & ys <= max(E_nodes(:,2)) & zs >= min(E_nodes(:,3)) & zs <= max(E_nodes(:,3));

        E_nodes  = local_nodes(ind,:);
        E_connect = [];
        % E_connect = delaunay(E_nodes(:,1), E_nodes(:,3));

        % % Renumber to global coordinates. This could be sped up
        % for iii = 1:size(E_connect,1)*size(E_connect,2)
        %     E_connect(iii) = ind(E_connect(iii));
        % end

    else % BELT ELECTRODE
        % Find all node indicies within a radius
        ind   = rangesearch(local_nodes, coord, E.E_rad);
        ind   = ind{1}';

        % Find all node coords based on indicies
        E_nodes = zeros(length(ind),3);
        for iii = 1:length(ind)
            E_nodes(iii,:) = local_nodes(ind(iii),:);
        end

        % % Create a trianagulation on the surface based on location of the nodes
        % [~,smallest_std] = min(std(E_nodes));
        % if smallest_std == 1
        %     E_connect = delaunay(E_nodes(:,2), E_nodes(:,3));
        % elseif smallest_std == 2
        %     E_connect = delaunay(E_nodes(:,1), E_nodes(:,3));
        % elseif smallest_std == 3
        %     E_connect = delaunay(E_nodes(:,1), E_nodes(:,2));
        % end
        % 
        % % Renumber to global coordinates. This could be sped up
        % for iii = 1:size(E_connect,1)*size(E_connect,2)
        %     E_connect(iii) = ind(E_connect(iii));
        % end

        E_connect = [];

        E_nodes = E_nodes;

    end
end

function [nodes, connect] = create_rows(local_nodes, coord, E)
    %{
    Create an row of electrodes and find the nodes/connectivity
    9/27/24 - Kyler Howard

    param: local_nodes - All nodes on to look at for making the electrode
    param: coord       - Reference coordinate for each electrode
                         Bottom left corner  for patches
                         Center of electrode for belt
    param: E           - Structure containing the size and type of electrode

    return: corners - Coordinates of the electrodes
    return: connect - Delauny trainagulation for the electrodes
    %}

    nodes   = cell(E.E_count(2),1);
    connect = cell(E.E_count(2),1);

    % Find the third electrode
    [nodes{3}, connect{3}] = create_electrode(local_nodes, coord, E);

    % Set up the distance between electrodes
    dist = E.E_width + E.gap_width;

    % Find the fourth electrode
    max_vals = max(nodes{3});
    coord(1) = max_vals(1) + E.gap_width/2;
    [nodes{4}, connect{4}] = create_electrode(local_nodes, coord, E);

    % Find the second electrode
    min_vals = min(nodes{3});
    coord(1) = min_vals(1) - dist - E.gap_width/2;
    [nodes{2}, connect{2}] = create_electrode(local_nodes, coord, E);

    % Find the first electrode
    min_vals = min(nodes{2});
    coord(1) = min_vals(1) - dist - E.gap_width/2;
    [nodes{1}, connect{1}] = create_electrode(local_nodes, coord, E);

    % % Find the other 3 electrodes
    % [corners(13:16, :), connect{4}] = create_electrode(local_nodes, coord +   [dist,0,0], E);
    % [corners(5:8, :),   connect{2}] = create_electrode(local_nodes, coord -   [dist,0,0], E);
    % [corners(1:4, :),   connect{1}] = create_electrode(local_nodes, coord - 2*[dist,0,0], E);
end

function [E_nodes, connect] = create_patch(local_nodes, coord, E, FoB)
    %{
    Create a patch of electrodes and find the nodes/connectivity
    9/27/24 - Kyler Howard

    param: local_nodes - All nodes on to look at for making the electrode
    param: coord       - Reference coordinate for each electrode
                         Bottom left corner  for patches
                         Center of electrode for belt
    param: E           - Structure containing the size and type of electrode
    param: FoB         - Flag for if we are looking at the front or back of the body

    return: corners - Coordinates of the electrodes
    return: connect - Delauny trainagulation for the electrodes
    %}

    E_nodes      = cell(E.E_count(1) * E.E_count(2), 1);
    row_nodes    = cell(E.E_count(1), 1);
    temp_connect = cell(E.E_count(1),1);

    % Find the bottom row
    [row_nodes{4}, temp_connect{4}] = create_rows(local_nodes, coord, E);
    for i = 1:4
        E_nodes{12+i} = row_nodes{4}{i};
    end

    % Set up the distance between rows
    % dist = E.E_height + E.gap_height;
    dist = E.E_height/2 + 4/3*E.gap_height;

    % Find the second row
    % max(corners(49:64,3))
    max_vals = max([max(E_nodes{13}); max(E_nodes{14}); max(E_nodes{15}); max(E_nodes{16})]);
    [coord2, ~]= find_center(local_nodes, max_vals(3) + dist, FoB);
    % coord2 = find_center(local_nodes, coord(3) + dist, FoB);
    [row_nodes{3},  temp_connect{3}] = create_rows(local_nodes, coord2, E);
    for i = 1:4
        E_nodes{8+i} = row_nodes{3}{i};
    end

    % Find the third row
    % max(corners(33:48,3))
    max_vals = max([max(E_nodes{9}); max(E_nodes{10}); max(E_nodes{11}); max(E_nodes{12})]);
    [coord3, ~] = find_center(local_nodes, max_vals(3) + dist, FoB);
    % coord3 = find_center(local_nodes, coord2(3) + dist, FoB);
    [row_nodes{2},  temp_connect{2}] = create_rows(local_nodes, coord3, E);
    for i = 1:4
        E_nodes{4+i} = row_nodes{2}{i};
    end

    % Find the fourth row
    % max(corners(17:32,3))
    max_vals = max([max(E_nodes{5}); max(E_nodes{6}); max(E_nodes{7}); max(E_nodes{8})]);
    [coord4, ~] = find_center(local_nodes, max_vals(3) + dist, FoB);
    % coord4 = find_center(local_nodes, coord3(3) + dist, FoB);
    [row_nodes{1},  temp_connect{1}] = create_rows(local_nodes, coord4, E);
    for i = 1:4
        E_nodes{i} = row_nodes{1}{i};
    end

    % Reformat the output
    connect = cell(16,1);
    for i = 1:4
        connect{i}    = temp_connect{1}{i};
        connect{i+4}  = temp_connect{2}{i};
        connect{i+8}  = temp_connect{3}{i};
        connect{i+12} = temp_connect{4}{i};
    end
end

function [E_nodes, E_connect] = create_belt(local_nodes, E_plane, E)
    %{
    Create a belt of electrodes and find the nodes/connectivity
    9/27/24 - Kyler Howard

    param: local_nodes - All nodes on to look at for making the electrode
    param: E_plane     - A subset of nodes to look through to find centers of electrodes quickly
    param: E           - Structure containing the size and type of electrode

    return: E_nodes   - Coordinates of the electrodes
    return: E_connect - Delauny trainagulation for the electrodes
    %}

    E_nodes   = cell(16, 1);
    E_connect = cell(16,1);

    center = (max(E_plane,[],1) + min(E_plane,[],1)) / 2;

    i = 1;
    for theta = 3*pi/2 : -(2*pi)/16 : -pi/2 + (2*pi)/16 % Shift to start at the back of the body
        % 3*pi/2 : -(2*pi)/16 : -pi/2 + (2*pi)/16
        radius =  Parameratize_Bdry(E_plane, 15, theta);
        c_point = [center(1) + radius*cos(theta), center(2) + radius*sin(theta), center(3)];

        E_center = E_plane(dsearchn(E_plane,c_point),:);
        
        [E_nodes{i}, E_connect{i}] = create_electrode(local_nodes, E_center, E);

        % figure
        % hold on
        % scatter3(nodes{i}(:,1),nodes{i}(:,2),nodes{i}(:,3),'r','filled')
        % trimesh(connect{i}, nodes{i}(:,1),nodes{i}(:,2),nodes{i}(:,3))

        i = i + 1;
    end

end


function [E_nodes, perim_mm_high] = Make_Electrodes3(boundary_nodes, all_nodes, body_faces,  sbj_info, flags)
    %{
    Find the nodes that make up all 32 electrodes, for either a belt or patch configuration
    Updated to compute area of each electrode, and make sure they are the correct sizes
    1/29/25 - Kyler Howard

    param: boundary_nodes - All nodes on the boundary of the surface
    param: sbj_info       - Heights of anatomical markers for electrode placement
    param: flags          - Various flags controlling plotting and other parameters

    return: E_nodes   - 1x32 Cell array containing electrode nodes
    return: perim_mm  - Double of the perimeter around the body in mm
    %}

% ----------------------------------------------------------------------- %
%% ------------------------------- Setup -------------------------------- %
% ----------------------------------------------------------------------- %
    % 1: Large patch front back
    % 2: Small patch front back
    % 3: Two rows of large belts
    % 4: Two rows of small belts

    % Constructing large GE patch
    L_square.type       = "patch";
    L_square.shape      = "rectangle";
    L_square.E_width    = 10;               % mm
    L_square.E_height   = 10;               % mm
    L_square.E_area     = L_square.E_width * L_square.E_height;
    L_square.E_count    = [4,4];            % Electrodes per row and per column
    L_square.gap_width  = 2.5;              % mm (edge to edge)
    L_square.gap_height = 2.5;              % mm (edge to edge)
    
    % Constructing small GE patch
    S_square.type       = "patch";
    S_square.shape      = "rectangle";
    S_square.E_width    = 7;                % mm
    S_square.E_height   = 7;                % mm
    S_square.E_area     = S_square.E_width * S_square.E_height;
    S_square.E_count    = [4,4];            % Electrodes per row and per column
    S_square.gap_width  = 2.5;              % mm (edge to edge)
    S_square.gap_height = 2.5;              % mm (edge to edge)
    
    % Constructing large GE belt
    L_belt.type    = "belt";
    L_belt.shape   = "circle";
    L_belt.E_dia   = 17;                    % mm
    L_belt.E_rad   = L_belt.E_dia / 2;      % mm
    L_belt.E_area  = pi * L_belt.E_rad^2;   % mm²
    L_belt.E_count = 16;                     % Electrodes per row
    
    % Constructing small GE belt
    S_belt.type    = "belt";
    S_belt.shape   = "circle";
    S_belt.E_dia   = 12;                    % mm
    S_belt.E_rad   = S_belt.E_dia / 2;      % mm
    S_belt.E_area  = pi * S_belt.E_rad^2;   % mm²
    S_belt.E_count = 16;                     % Electrodes per row

    % Constructing custom electrode setup
    E_custom.type  = flags.E_type;
    E_custom.shape = flags.E_shape;
    if E_custom.type == "patch"
        E_custom.E_count = flags.E_count;
        E_custom.gap_width = flags.gap_width;
        E_custom.gap_height = flags.gap_height;
    elseif E_custom.type == "belt"
        E_custom.E_count = flags.E_count;
    end

    if E_custom.shape == "circle"
        E_custom.E_dia  = flags.E_dia;
        E_custom.E_rad  = flags.E_dia / 2;
        E_custom.E_area = pi * E_custom.E_rad^2;
    elseif E_custom.shape == "rectangle"
        E_custom.E_width  = flags.E_width;
        E_custom.E_height = flags.E_height;
        E_custom.E_area   = E_custom.E_width * E_custom.E_height;
    end

    choices = {L_square, S_square, L_belt, S_belt, E_custom};
    E = choices{flags.E_choice};
    
% ----------------------------------------------------------------------- %
%% --------------------------- Center Nodes ----------------------------- %
% ----------------------------------------------------------------------- %
    % Find the plane in which the center lies
    if E.type == "patch"
        try
            plane_gap = E.E_rad;
        catch
            plane_gap = E.E_height / 2;
        end

        z_check = (boundary_nodes(:,3) > sbj_info.carina - plane_gap) & (boundary_nodes(:,3) < sbj_info.carina + plane_gap);
        plane_high = boundary_nodes(z_check,:);

        % Find the nodes that are the center of the front/back of the body
        x_mid = (min(plane_high(:,1)) + max(plane_high(:,1))) / 2;
        y_min = min(plane_high(:,2));
        y_max = max(plane_high(:,2));
        front = plane_high(dsearchn(plane_high, [x_mid, y_max, sbj_info.carina]),:);
        back  = plane_high(dsearchn(plane_high, [x_mid, y_min, sbj_info.carina]),:);
        % front = find_node(plane, x_mid, y_min, sbj_info.carina);
        % back  = find_node(plane, x_mid, y_max, sbj_info.carina);

    elseif E.type == "belt"
        if E.shape == "circle"
            plane_gap = E.E_rad;
        elseif E.shape == "rectangle"
            plane_gap = E.E_height;
        end

        z_check_high = (boundary_nodes(:,3) > sbj_info.T5 - plane_gap) & (boundary_nodes(:,3) < sbj_info.T5 + plane_gap);
        plane_high   = boundary_nodes(z_check_high,:);

        % height_low  = sbj_info.carina - E.E_dia - E.gap;
        z_check_low = (boundary_nodes(:,3) > sbj_info.T8 - plane_gap) & (boundary_nodes(:,3) < sbj_info.T8 + plane_gap);
        plane_low   = boundary_nodes(z_check_low,:);
    end

    

% ----------------------------------------------------------------------- %
%% ------------------------- Finding Perimeter -------------------------- %
% ----------------------------------------------------------------------- %
    center_high = (max(plane_high,[],1) + min(plane_high,[],1)) / 2;
    if E.type == "belt"
        center_low = (max(plane_low,[],1) + min(plane_low,[],1)) / 2;
    end


    i = 1;
    point_low     = zeros(200, 3);
    point_high    = zeros(200, 3);
    perim_mm_low  = 0;
    perim_mm_high = 0;
    for theta = 0 : (2*pi)/200 : 2*pi - (2*pi)/200 
        radius_high =  Parameratize_Bdry(plane_high, 15, theta);

        point_high(i,:) = [center_high(1) + radius_high*cos(theta), center_high(2) + radius_high*sin(theta), center_high(3)];

        if i >= 2  
            dist_high = sqrt((point_high(i-1,1) - point_high(i,1))^2 + (point_high(i-1,2) - point_high(i,2))^2);
            perim_mm_high = perim_mm_high + dist_high;
        end

        % Repeat for the lower belt
        if E.type == "belt"
            radius_low =  Parameratize_Bdry(plane_low, 15, theta);

            point_low(i,:) = [center_low(1) + radius_low*cos(theta), center_low(2) + radius_low*sin(theta), center_low(3)];
            
            if i >= 2  
                dist_low = sqrt((point_low(i-1,1) - point_low(i,1))^2 + (point_low(i-1,2) - point_low(i,2))^2);
                perim_mm_low = perim_mm_low + dist_low;
            end
        end

        i = i + 1;
    end

    if flags.plot_slices && flags.plot_electrodes
        figure()
            hold on
            scatter(plane_high(:,1), plane_high(:,2))
            plot(point_low(:,1),   point_low(:,2), 'r', 'linewidth', 1.5)
            legend("Exact Points", "Parameratized Boundary", 'location','southoutside')
    end

% ----------------------------------------------------------------------- %
%% ---------------------- Finding Electrode Nodes ----------------------- %
% ----------------------------------------------------------------------- %
    if E.type == "patch"
        % Find the corners of the electrodes
        front_nodes = create_patch(boundary_nodes, front, E, "front", all_nodes, body_faces, flags);
        back_nodes  = create_patch(boundary_nodes, back,  E, "back", all_nodes, body_faces, flags);
        E_nodes     = [front_nodes; back_nodes];

    elseif E.type == "belt"
        top_nodes = create_belt(boundary_nodes, plane_high, E, all_nodes, body_faces, flags, perim_mm_high);
        bot_nodes = create_belt(boundary_nodes, plane_low,  E, all_nodes, body_faces, flags, perim_mm_low);
        E_nodes   = cat(1, bot_nodes,   top_nodes)';
    end

    % Reorder electrode placement if using the 4x8 pattern
    if flags.CP_choice == 2
        fprintf("   Placing 4x8 Electrode Arrays\n")
        all_ind = [13,14,15,16,32,31,30,29,9,10,11,12,28,27,26,25,5,6,7,8,24,23,22,21,1,2,3,4,20,19,18,17];
        E_nodes = E_nodes(all_ind);
    end


% ----------------------------------------------------------------------- %
%% ------------------------------ Plotting ------------------------------ %
% ----------------------------------------------------------------------- %

    if flags.plot_electrodes == 1 && E.type == "patch"
        figure(); 
            hold on; 
            for cell_i = 1:size(E_nodes,1) / 2
                scatter3(E_nodes{cell_i}(:,1), E_nodes{cell_i}(:,2), E_nodes{cell_i}(:,3), 'filled', 'o');
            end
            for cell_i = (size(E_nodes,1) / 2) + 1:size(E_nodes,1)
                scatter3(E_nodes{cell_i}(:,1), E_nodes{cell_i}(:,2), E_nodes{cell_i}(:,3), 'filled', 'square');
            end
            legend('Location','eastoutside')
            xlabel('X (mm)')
            ylabel('Y (mm)')
            zlabel('Z (mm)')
            axis equal
    elseif flags.plot_electrodes == 1 && E.type == "belt"
        figure(); 
            hold on
                scatter3(E_nodes{1}(:,1), E_nodes{1}(:,2), E_nodes{1}(:,3), 'filled', 'square')
            for node_i = 2:length(E_nodes)
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

function E_nodes = create_electrode(local_nodes, coord, E, all_nodes, body_faces, flags)
    %{
    Create an individual electrode and find the nodes/connectivity
    9/27/24 - Kyler Howard

    param: local_nodes - All nodes on to look at for making the electrode
    param: coord       - Reference coordinate for each electrode, always the center
    param: E           - Structure containing the size and type of electrode

    return: corners   - Coordinates of the electrode
    %}

    good_electrode = 0;
    search_dist    = 0;
    while good_electrode == 0
        % Find all node indicies within a radius
        if E.shape == "circle"
            ind = rangesearch(local_nodes, coord, E.E_rad + search_dist);
            ind = ind{1}';
        elseif E.shape == "rectangle"
            % Until I can make the electrode touching, go half left/right/up
            xcheck = (local_nodes(:,1) >= coord(1) - (E.E_width + search_dist/2)/2)  & (local_nodes(:,1) <= coord(1) + (E.E_width + search_dist/2)/2);
            if E.type == "belt"
            ycheck = (local_nodes(:,2) >= coord(2) - (E.E_width + search_dist/2)/2)  & (local_nodes(:,2) <= coord(2) + (E.E_width + search_dist/2)/2);
            elseif E.type == "patch"
                ycheck = (local_nodes(:,2) >= coord(2) - (E.E_width + search_dist/2))  & (local_nodes(:,2) <= coord(2) + (E.E_width + search_dist/2));
            end
            % zcheck = (local_nodes(:,3) >= coord(3) - (E.E_height + 2*search_dist)/2) & (local_nodes(:,3) <= coord(3) + (E.E_height + 2*search_dist)/2);
            zcheck = (local_nodes(:,3) >= coord(3) - (E.E_height + search_dist)/2) & (local_nodes(:,3) <= coord(3) + (E.E_height + search_dist)/2);
            binary_check = xcheck & ycheck & zcheck;

            ind = zeros(sum(binary_check), 1);
            j = 1;
            for i = 1:size(binary_check,1)
                if binary_check(i) == 1
                    ind(j) = i;
                    j      = j + 1;
                end
            end
        end

        % Find all node coords based on indicies
        E_nodes = zeros(length(ind),3);
        for iii = 1:length(ind)
            E_nodes(iii,:) = local_nodes(ind(iii),:);
        end

        E_connect = Align_Electrode_Faces2(all_nodes, body_faces, E_nodes);

        E_area = 0;
        for j = 1:size(E_connect,1)
            face_nodes = E_connect(j,:);
            point1 = all_nodes(face_nodes(1),:);
            point2 = all_nodes(face_nodes(2),:);
            point3 = all_nodes(face_nodes(3),:);
            E_area = E_area + norm(cross(point3-point1, point3-point2))/2;
        end

        if E_area < E.E_area
            good_electrode = 0;
            search_dist = search_dist + 0.5;
        else
            good_electrode = 1;
        end

    end
end

function nodes = create_rows(local_nodes, coord, E, all_nodes, body_faces, flags)
    %{
    Create an row of electrodes and find the nodes/connectivity
    9/27/24 - Kyler Howard

    param: local_nodes - All nodes on to look at for making the electrode
    param: coord       - Reference coordinate for each electrode
                         Bottom left corner  for patches
                         Center of electrode for belt
    param: E           - Structure containing the size and type of electrode

    return: corners - Coordinates of the electrodes
    %}

    nodes   = cell(E.E_count(2),1);

    try
        coord = coord + [E.E_rad + E.gap_width/2, 0, 0];
    catch
        coord = coord + [E.E_width/2 + E.gap_width/2, 0, 0];
    end

    % Find the third electrode
    nodes{3} = create_electrode(local_nodes, coord, E, all_nodes, body_faces, flags);

    % Set up the distance between electrodes
    try
        dist = E.E_rad + E.gap_width;
        % dist = E.E_rad + 5/3*E.gap_width;
    catch
        % dist = E.E_width/2 + E.gap_width;
        dist = E.E_width/2 + 5/3*E.gap_width;
    end

    % Find the fourth electrode
    max_vals = max(nodes{3});
    coord(1) = max_vals(1) + dist;
    nodes{4} = create_electrode(local_nodes, coord, E, all_nodes, body_faces, flags);

    % Find the second electrode
    min_vals = min(nodes{3});
    coord(1) = min_vals(1) - dist;
    nodes{2} = create_electrode(local_nodes, coord, E, all_nodes, body_faces, flags);


    % Find the first electrode
    min_vals = min(nodes{2});
    coord(1) = min_vals(1) - dist;
    nodes{1} = create_electrode(local_nodes, coord, E, all_nodes, body_faces, flags);

    % % Find the other 3 electrodes
    % corners(13:16, :) = create_electrode(local_nodes, coord +   [dist,0,0], E);
    % corners(5:8, :) = create_electrode(local_nodes, coord -   [dist,0,0], E);
    % corners(1:4, :) = create_electrode(local_nodes, coord - 2*[dist,0,0], E);
end

function E_nodes = create_patch(local_nodes, start_point, E, FoB, all_nodes, body_faces, flags)
    %{
    Create a patch of electrodes and find the nodes/connectivity
    9/27/24 - Kyler Howard

    param: local_nodes - All nodes on to look at for making the electrode
    param: coord       - Reference coordinate for each electrode
                         Bottom left corner  for patches
    param: E           - Structure containing the size and type of electrode
    param: FoB         - Flag for if we are looking at the front or back of the body

    return: corners - Coordinates of the electrodes
    %}

    E_nodes   = cell(E.E_count(1) * E.E_count(2), 1);
    row_nodes = cell(E.E_count(1), 1);

    % Find the third row
    row_nodes{3} = create_rows(local_nodes, start_point, E, all_nodes, body_faces, flags);
    for i = 1:4
        E_nodes{8+i} = row_nodes{3}{i};
    end

    % Set up the distance between rows
    try
        dist = E.E_rad + E.gap_height;
        % dist = E.E_rad + 5/3*E.gap_height;
    catch
        % dist = E.E_height/2 + E.gap_height;
        dist = E.E_height/2 + 5/3*E.gap_height;
    end

    % Find the fourth row
    % max(corners(49:64,3))
    min_vals = min([min(E_nodes{9}); min(E_nodes{10}); min(E_nodes{11}); min(E_nodes{12})]);
    [coord2, ~]= find_center(local_nodes, min_vals(3) - dist, FoB);
    % coord2 = find_center(local_nodes, coord(3) + dist, FoB);
    row_nodes{4} = create_rows(local_nodes, coord2, E, all_nodes, body_faces, flags);
    for i = 1:4
        E_nodes{12+i} = row_nodes{4}{i};
    end

    % Find the second row
    % max(corners(33:48,3))
    max_vals = max([max(E_nodes{9}); max(E_nodes{10}); max(E_nodes{11}); max(E_nodes{12})]);
    [coord3, ~] = find_center(local_nodes, max_vals(3) + dist, FoB);
    % coord3 = find_center(local_nodes, coord2(3) + dist, FoB);
    row_nodes{2} = create_rows(local_nodes, coord3, E, all_nodes, body_faces, flags);
    for i = 1:4
        E_nodes{4+i} = row_nodes{2}{i};
    end

    % Find the first row
    % max(corners(17:32,3))
    max_vals = max([max(E_nodes{5}); max(E_nodes{6}); max(E_nodes{7}); max(E_nodes{8})]);
    [coord4, ~] = find_center(local_nodes, max_vals(3) + dist, FoB);
    % coord4 = find_center(local_nodes, coord3(3) + dist, FoB);
    row_nodes{1} = create_rows(local_nodes, coord4, E, all_nodes, body_faces, flags);
    for i = 1:4
        E_nodes{i} = row_nodes{1}{i};
    end
end

function E_nodes = create_belt(local_nodes, E_plane, E, all_nodes, body_faces, flags, perim_mm)
    %{
    Create a belt of electrodes and find the nodes/connectivity
    9/27/24 - Kyler Howard

    param: local_nodes - All nodes on to look at for making the electrode
    param: E_plane     - A subset of nodes to look through to find centers of electrodes quickly
    param: E           - Structure containing the size and type of electrode

    return: E_nodes   - Coordinates of the electrodes
    %}

    E_nodes   = cell(E.E_count, 1);

    center = (max(E_plane,[],1) + min(E_plane,[],1)) / 2;

    % KH: Equal Angle Electrodes
    % i = 1;
    % for theta = pi/2 : -(2*pi)/E.E_count : -(3*pi)/2 + (2*pi)/E.E_count % Shift to start at the back of the body
    %     % 3*pi/2 : -(2*pi)/16 : -pi/2 + (2*pi)/16
    %     % Parameratize using 15 terms of a Fourier expansion
    %     radius =  Parameratize_Bdry(E_plane, 15, theta);
    %     c_point = [center(1) + radius*cos(theta), center(2) + radius*sin(theta), center(3)];
    % 
    %     E_center = E_plane(dsearchn(E_plane,c_point),:);
    % 
    %     E_nodes{i} = create_electrode(local_nodes, E_center, E, all_nodes, body_faces, flags);
    % 
    %     i = i + 1;
    % end

    % KH: Equal Arc Length Electrodes
    i = 1;
    j = 1;
    arc_length = 0;
    point = zeros(200, 3);
    for theta = 3*pi/2 : -(2*pi)/200 : -pi/2 + (2*pi)/200
        radius =  Parameratize_Bdry(E_plane, 15, theta);
        point(j,:) = [center(1) + radius*cos(theta), center(2) + radius*sin(theta), center(3)];
        
        % Always make electrode on the back
        if theta == 3*pi/2
            c_point    = point(j,:);
            E_center   = E_plane(dsearchn(E_plane,c_point),:);
            E_nodes{i} = create_electrode(local_nodes, E_center, E, all_nodes, body_faces, flags);
            i = i + 1;
        end

        % Measure the distance between sweeps
        if j >= 2  
            dist = sqrt((point(j-1,1) - point(j,1))^2 + (point(j-1,2) - point(j,2))^2);
            arc_length = arc_length + dist;
        end

        % Check if we have moved around enough
        % if arc_length >= perim_mm / (E.E_count + 1)
        % try
        %     goal_arc_length = (perim_mm - 2*E.E_width) / E.E_count;
        % catch
            % goal_arc_length = (perim_mm - 2*E.E_rad) / E.E_count;
            goal_arc_length = perim_mm / (E.E_count + 0.7);
        % end

        if arc_length >= goal_arc_length
            % Reset arc length
            arc_length = 0;

            % Find the electrode
            c_point    = point(j,:);
            E_center   = E_plane(dsearchn(E_plane,c_point),:);
            E_nodes{i} = create_electrode(local_nodes, E_center, E, all_nodes, body_faces, flags);
            i = i + 1;
        end

        % Update the point index
        j = j + 1;

        if i == 17
            break
        end
    end

end


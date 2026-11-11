function [E_nodes, perim_mm, vert_gap] = Make_Electrodes(boundary_nodes, all_nodes, body_faces, sbj_info, flags)
    %{
    Find the nodes that make up all L electrodes, for either a belt or patch configuration
    Updated to compute area of each electrode, and make sure they are the correct sizes
    1/29/25 - Kyler Howard

    param: boundary_nodes - All nodes on the boundary of the surface
    param: all_nodes      - All nodes in the entire body
    param: body_faces     - Triangular faces that make up the surface of the body
    param: sbj_info       - Heights of anatomical markers for electrode placement
    param: flags          - Various flags controlling plotting and other parameters

    return: E_nodes       - 1xL Cell array containing electrode nodes
    return: perim_mm_high - The perimeter around the body in mm
    %}

% ----------------------------------------------------------------------- %
%% ------------------------------- Setup -------------------------------- %
% ----------------------------------------------------------------------- %
    % Pull the most recent electrode settings
    [E, ~] = Construct_Electrode_Settings(flags);
    
% ----------------------------------------------------------------------- %
%% --------------------------- Center Nodes ----------------------------- %
% ----------------------------------------------------------------------- %
    % Find the plane in which the center lies
    if E.type == "patch"
        nheights            = 1;
        heights             = sbj_info.carina;
        [front, E_plane{1}] = Find_Center(boundary_nodes, heights, E, "front");
        [back,  ~]          = Find_Center(boundary_nodes, heights, E, "back");
    elseif E.type == "belt"
        nheights = size(E.E_count,2);
        E_plane = cell(1,nheights);
        if nheights == 1
            heights         = mean([sbj_info.T5, sbj_info.T8]);
            [~, E_plane{1}] = Find_Center(boundary_nodes,  heights, E);
        elseif nheights == 2
            heights         = [sbj_info.T5, sbj_info.T8];
            [~, E_plane{1}] = Find_Center(boundary_nodes, heights(1), E);
            [~, E_plane{2}] = Find_Center(boundary_nodes, heights(2), E);
        else
            error("%d rows of electrodes are unsupported\n", nheights)
        end
    end

% ----------------------------------------------------------------------- %
%% ------------------------- Finding Perimeter -------------------------- %
% ----------------------------------------------------------------------- %
    if flags.verbose
        make_param_string = '      Parameterizing the boundary...';
        fprintf("%s", make_param_string)
    end
    param_start   = tic;
    body_center   = cell(1,nheights);
    bdry_parm_pts = cell(1,nheights);
    perim_mm      = cell(1,nheights);
    for i = 1:nheights
        % Find the center of the body for each electrode height
        body_center{i} = (max(E_plane{i},[],1) + min(E_plane{i},[],1)) / 2;

        % Parameratize
        j = 1;
        n_points         = 2000;
        param_terms      = 40;
        bdry_parm_pts{i} = zeros(n_points, 3);
        perim_mm{i}      = 0;
        for theta = 0 : (2*pi)/n_points : 2*pi - (2*pi)/n_points 
            % Find the radius based on Fourier parameterization at set theta
            radius =  Parameratize_Bdry(E_plane{i}, param_terms, theta);

            % Calculate parameratized points
            bdry_parm_pts{i}(j,:) = [body_center{i}(1) + radius*cos(theta), body_center{i}(2) + radius*sin(theta), body_center{i}(3)];

            % If this isn't the first point, calculate the distance traveled between points
            if j >= 2  
                dist        = norm(bdry_parm_pts{i}(j-1,:) - bdry_parm_pts{i}(j,:));
                perim_mm{i} = perim_mm{i}  + dist;
            end

            j = j + 1;
        end
    end
    param_stop   = toc(param_start);

    if flags.verbose == 1
        for ii = 1:length(make_param_string)
            fprintf("\b")
        end
        fprintf("      Parameterized the boundary in %.2f seconds\n", param_stop)
    end

% ----------------------------------------------------------------------- %
%% ---------------------- Finding Electrode Nodes ----------------------- %
% ----------------------------------------------------------------------- %
    if flags.verbose
        make_elec_string = '      Making electrodes...';
        fprintf("%s", make_elec_string)
    end
    if E.type == "patch"
        elec_start = tic;
        % Adjust boundary nodes to not include the top/bottom plane. Remove top/bot 0.5 mm
        tube_nodes = boundary_nodes(boundary_nodes(:,3)>0.5 & boundary_nodes(:,3)<max(boundary_nodes(:,3))-0.5, :);
        
        % Only look at the nodes near the patch. Don't waste time looking at nodes near the armpits
        front_tube_nodes   = tube_nodes(tube_nodes(:,2) > mean(tube_nodes(:,2)),:);
        back_tube_nodes    = tube_nodes(tube_nodes(:,2) < mean(tube_nodes(:,2)),:);
        if mod(E.E_count(2), 2) == 0
            gap = (E.E_count(2)/2 + 1)*E.E_width + (E.E_count(2)/2)*E.gap_width;
        else
            gap = ceil(E.E_count(2)/2)*E.E_width + ceil(E.E_count(2)/2)*E.gap_width;
        end
        middle_front_nodes = front_tube_nodes(front_tube_nodes(:,1) > mean(front_tube_nodes(:,1)) - gap & front_tube_nodes(:,1) < mean(front_tube_nodes(:,1)) + gap,:);
        middle_back_nodes  = back_tube_nodes(back_tube_nodes(:,1) > mean(back_tube_nodes(:,1)) - gap & back_tube_nodes(:,1) < mean(back_tube_nodes(:,1)) + gap,:);

        % Find the corners of the electrodes
        front_nodes = create_patch(middle_front_nodes, front, E, "front", all_nodes, body_faces);
        back_nodes  = create_patch(middle_back_nodes,  back,  E, "back",  all_nodes, body_faces);
        E_nodes     = [front_nodes; back_nodes];

        % Find vertical gap
        n_gaps = flags.E_count(1) - 1;
        vert_gap = zeros(n_gaps, 1);
        for ii = 1:n_gaps
            nEl_p_row = flags.E_count(2);
            % Extract top nodes from front/back
            bot_nodes = [front_nodes(nEl_p_row*(ii-1) + 1 : nEl_p_row*ii); back_nodes(nEl_p_row*(ii-1) + 1 : nEl_p_row*ii)];
            top_nodes = [front_nodes(nEl_p_row*ii + 1 : nEl_p_row*(ii+1)); back_nodes(nEl_p_row*ii + 1 : nEl_p_row*(ii+1))];

            meanTop = mean(cellfun(@(A) mean(A(:,3)), top_nodes));
            meanBot = mean(cellfun(@(A) mean(A(:,3)), bot_nodes));
            vert_gap(ii) = abs(meanTop - meanBot)/10; % Center-to-center gap in cm

            meanTopFront = mean(cellfun(@(A) mean(A(:,3)), front_nodes(nEl_p_row*ii + 1 : nEl_p_row*(ii+1))));
            meanTopBack  = mean(cellfun(@(A) mean(A(:,3)), back_nodes(nEl_p_row*ii + 1 : nEl_p_row*(ii+1))));
            meanBotFront = mean(cellfun(@(A) mean(A(:,3)), front_nodes(nEl_p_row*(ii-1) + 1 : nEl_p_row*ii)));
            meanBotBack  = mean(cellfun(@(A) mean(A(:,3)), back_nodes(nEl_p_row*(ii-1) + 1 : nEl_p_row*ii)));
            vert_gap(ii) = mean([abs(meanTopFront - meanBotFront)/10, abs(meanTopBack - meanBotBack)/10]); % Center-to-center gap in cm
        end
        vert_gap = mean(vert_gap);
        elec_stop = toc(elec_start);

    elseif E.type == "belt"
        elec_start = tic;
        temp_E_nodes = create_belts(boundary_nodes, E_plane, E, all_nodes, body_faces, flags, perim_mm);
        % bot_nodes = create_belts(boundary_nodes, plane_low,  E, all_nodes, body_faces, flags, perim_mm_low);
        % E_nodes   = cat(1, bot_nodes,   top_nodes)';
        E_nodes = [];
        for i = 1:nheights
            E_nodes = cat(1, E_nodes, temp_E_nodes{i});
        end
        E_nodes = E_nodes';

        % Find vertical gap
        meanZs = zeros(1,nheights);
        for i = 1:nheights
            meanZs(i) = mean(cellfun(@(A) mean(A(:,3)), temp_E_nodes{i}));
        end

        % Find the average of adjacent absolute differences of row heights in cm
        % Center-to-center gap
        if nheights == 1
            vert_gap = meanZs/10;
        else
            abs_difs = abs(diff(meanZs));
            vert_gap = mean(abs_difs)/10;
        end
        elec_stop = toc(elec_start);
    end
    if flags.verbose
        for ii = 1:length(make_elec_string)
            fprintf("\b")
        end
        fprintf("      Made electrodes in %.2f seconds\n", elec_stop)
    end

    % Reorder electrode placement if using the 4x8 pattern
    if flags.CP_choice == 2 && E.type == "patch"
        fprintf("   Placing 4x8 electrode arrays\n")
        % 1/1/26 KH: GE updated the order for the patch array
        % new_ind = [13,14,15,16,32,31,30,29,9,10,11,12,28,27,26,25,5,6,7,8,24,23,22,21,1,2,3,4,20,19,18,17]; % KH 12/5/25 - GE Changed the 4x8 patch configuration
        new_ind = [31,32,16,15,14,13,29,30,27,28,12,11,10,9,25,26,23,24,8,7,6,5,21,22,19,20,4,3,2,1,17,18];
        E_nodes = E_nodes(new_ind);
    end

    % Find average perimeter
    perim_mm = mean(cat(3, perim_mm{:}));
% ----------------------------------------------------------------------- %
%% ------------------------------ Plotting ------------------------------ %
% ----------------------------------------------------------------------- %

    if flags.plot_electrodes == 1 && E.type == "patch"
        figure(); 
            hold on; 
            for cell_i = 1:size(E_nodes,1) / 2
                if cell_i == 1
                    scatter3(E_nodes{cell_i}(:,1), E_nodes{cell_i}(:,2), E_nodes{cell_i}(:,3), 'k', 'filled', 'o');
                else
                    scatter3(E_nodes{cell_i}(:,1), E_nodes{cell_i}(:,2), E_nodes{cell_i}(:,3), 'filled', 'o');
                end
            end
            for cell_i = (size(E_nodes,1) / 2) + 1:size(E_nodes,1)
                if cell_i == (size(E_nodes,1) / 2) + 1
                    scatter3(E_nodes{cell_i}(:,1), E_nodes{cell_i}(:,2), E_nodes{cell_i}(:,3), 'k', 'filled', 'square');
                else
                    scatter3(E_nodes{cell_i}(:,1), E_nodes{cell_i}(:,2), E_nodes{cell_i}(:,3), 'filled', 'square');
                end
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
                if node_i == 1 || node_i == 17
                    scatter3(E_nodes{node_i}(:,1), E_nodes{node_i}(:,2), E_nodes{node_i}(:,3), 'k', 'filled', 'square')
                else
                    scatter3(E_nodes{node_i}(:,1), E_nodes{node_i}(:,2), E_nodes{node_i}(:,3), 'filled')
                end
            end
            legend('Location','eastoutside')
            xlabel('X (mm)')
            ylabel('Y (mm)')
            zlabel('Z (mm)')
            axis equal
            for i = 1:nheights
                plot3(bdry_parm_pts{i}(:,1),bdry_parm_pts{i}(:,2),bdry_parm_pts{i}(:,3), 'r')
            end
    end
end
% ----------------------------------------------------------------------- %
%% -------------------------- Custom Functions -------------------------- %
% ----------------------------------------------------------------------- %

function E_nodes = create_electrode(local_nodes, coord, E, all_nodes, body_faces)
    %{
    Create an individual electrode and find the nodes/connectivity
    9/27/24 - Kyler Howard

    param: local_nodes - Subset of nodes to look at for making the electrode
    param: coord       - Reference coordinate for each electrode, always the center
    param: E           - Structure containing the size and type of electrode
    param: all_nodes   - Matrix of all nodes to calculate electrode areas
    param: body_faces  - Matrix of surface faces to calculate electrode areas

    return: E_nodes    - Coordinates of the electrode
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
                % Allow all nodes in the y direction
                ycheck = ones(length(local_nodes),1);
                % ycheck = (local_nodes(:,2) >= coord(2) - (E.E_width + search_dist/2))  & (local_nodes(:,2) <= coord(2) + (E.E_width + search_dist/2));
            end
            zcheck = (local_nodes(:,3) >= coord(3) - (E.E_height + search_dist)/2) & (local_nodes(:,3) <= coord(3) + (E.E_height + search_dist)/2);
            % zcheck = (local_nodes(:,3) >= coord(3) - (E.E_height + search_dist/5)/2) & (local_nodes(:,3) <= coord(3) + (E.E_height + search_dist/5)/2);
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

        E_connect = Align_Individual_Electrode_Faces(all_nodes, body_faces, E_nodes);

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
            search_dist = search_dist + 0.25; % Look 0.2 mm further
        else
            good_electrode = 1;
        end

        % Create an escape condition for each electrode
        if search_dist > 25
            good_electrode = 1;
        end
    end
end

function E_nodes = create_rows(local_nodes, mid_col_center, E, all_nodes, body_faces)
    %{
    Create an row of electrodes and find the nodes/connectivity
    Electrode 1 is the largest x value, and electrode 4 is the smallest x value
    9/27/24 - Kyler Howard

    param: local_nodes - All nodes on to look at for making the electrode
    param: coord       - Reference coordinate for each electrode (center)
    param: E           - Structure containing the size and type of electrode

    return: E_nodes - all nodes within the electrode
    %}

    E_nodes   = cell(E.E_count(2),1);

    % Find the "middle" column that we should be appying electrodes onto. This is either 
    % the true middle for odd number of columns or the right of the two center columns.
    middle_col = floor(mean(1:E.E_count(2)));
    if mod(E.E_count(2), 2) == 0
        if isnan(E.E_width)
            mid_col_center(:,1) = mid_col_center(:,1) + E.E_rad + E.gap_width/2;
        elseif isnan(E.E_rad)
            mid_col_center(:,1) = mid_col_center(:,1) + E.E_width/2 + E.gap_width/2;
        end
    end

    % Set up the distance between rows
    if isnan(E.E_width)
        dist = E.E_rad + E.E_width;
    elseif isnan(E.E_rad)
        dist = E.E_width/2 + E.E_width;
    end

    % Find the middle electrode
    E_nodes{middle_col} = create_electrode(local_nodes, mid_col_center, E, all_nodes, body_faces);

    % Build the electrode on -x side of this
    center_coord = mid_col_center;
    for col_i = middle_col + 1 : E.E_count(2)
        % Space the electrode out and find it
        min_vals        = min(E_nodes{col_i - 1});                                               % find the min x of the previous electrode
        center_coord(1) = min_vals(1) - dist;                                            % shift the point over
        E_nodes{col_i}  = create_electrode(local_nodes, center_coord, E, all_nodes, body_faces); % make the electrode
        
        % Determine if adjacent electrodes are touching each other side-to-side
        search_dist = 0;
        while ~isempty(intersect(E_nodes{col_i}, E_nodes{col_i-1}, "rows"))
            search_dist     = search_dist + 1; % Look 1 mm further away
            center_coord(1) = min_vals(1) - dist - search_dist;
            E_nodes{col_i}  = create_electrode(local_nodes, center_coord, E, all_nodes, body_faces);
        end
    end

    % Build the electrode on +x side of this
    center_coord = mid_col_center;
    for col_i = middle_col - 1 : -1 : 1
        % Space the electrode out and find it
        max_vals        = max(E_nodes{col_i + 1});                                               % find the max x of the previous electrode
        center_coord(1) = max_vals(1) + dist;                                            % shift the point over
        E_nodes{col_i}  = create_electrode(local_nodes, center_coord, E, all_nodes, body_faces); % make the electrode
        
        % Determine if adjacent electrodes are touching each other side-to-side
        search_dist = 0;
        while ~isempty(intersect(E_nodes{col_i}, E_nodes{col_i+1}, "rows"))
            search_dist     = search_dist + 1; % Look 1 mm further away
            center_coord(1) = max_vals(1) + dist + search_dist;
            E_nodes{col_i}  = create_electrode(local_nodes, center_coord, E, all_nodes, body_faces);
        end
    end
end

function E_nodes = create_patch(local_nodes, mid_row_center, E, FoB, all_nodes, body_faces)
    %{
    Create a patch of electrodes and find the nodes/connectivity
    Row 1 is the largest z value, and row 4 is the smallest z value
    Rows are centered on carina
    9/27/24 - Kyler Howard

    param: local_nodes    - All nodes on to look at for making the electrode
    param: mid_row_center - Reference coordinate for the third row (center)
    param: E              - Structure containing the size and type of electrode
    param: FoB            - Flag for if we are looking at the front or back of the body

    return: corners - Coordinates of the electrodes
    %}

    E_nodes   = cell(E.E_count(1) * E.E_count(2), 1);
    row_nodes = cell(E.E_count(1), 1);

    % Find the "middle" row that we should be appying electrodes onto. This is either 
    % the true middle for odd number of rows or the top of the two center rows.
    middle_row = floor(mean(1:E.E_count(1)));
    if mod(E.E_count(1), 2) == 0
        if isnan(E.E_height)
            mid_row_center(:,3) = mid_row_center(:,3) + E.E_rad + E.gap_height/2;
        elseif isnan(E.E_rad)
            mid_row_center(:,3) = mid_row_center(:,3) + E.E_height/2 + E.gap_height/2;
        end
    end

    % Set up the distance between rows
    if isnan(E.E_height)
        dist = E.E_rad + E.gap_height;
    elseif isnan(E.E_rad)
        dist = E.E_height/2 + E.gap_height;
    end

    % Find the middle row
    row_nodes{middle_row} = create_rows(local_nodes, mid_row_center, E, all_nodes, body_faces);
    for i = 1:E.E_count(2)
        E_nodes{E.E_count(2)*(middle_row-1)+i} = row_nodes{middle_row}{i};
    end

    % Build the rows below this
    for row_i = middle_row + 1 : E.E_count(1)
        % Space the row out and find it
        min_vals         = min(cell2mat(cellfun(@(x) min(x), row_nodes{row_i - 1}, 'UniformOutput', false))); % find the bottom of the previous row
        [row_center, ~]  = Find_Center(local_nodes, min_vals(3) - dist, E, FoB);                              % find the closest center point after walking down the vert gap and half electrode height
        row_nodes{row_i} = create_rows(local_nodes, row_center, E, all_nodes, body_faces);                    % make the row
        
        % Set aside the electrode and check if rows are touching top-to-bottom
        touch_check = cell(1,E.E_count(2));
        for i = 1:E.E_count(2)
            % Save the individual electrode
            E_nodes{E.E_count(2)*(row_i-1)+i} = row_nodes{row_i}{i};

            % Check if they're touching
            touch_check{i} = intersect(row_nodes{row_i}{i}, row_nodes{row_i-1}{i},"rows");
        end

        % If touching, move further down
        search_dist = 0;
        while any(cellfun(@(x) ~isempty(x), touch_check))
            search_dist      = search_dist + 1; % Look 1 mm further away
            [row_center, ~]  = Find_Center(local_nodes, min_vals(3) - dist - search_dist, E, FoB); % find the closest center point after walking down the vert gap and half electrode height
            row_nodes{row_i} = create_rows(local_nodes, row_center, E, all_nodes, body_faces);     % make the row
            for i = 1:E.E_count(2)
                % Save the individual electrode
                E_nodes{E.E_count(2)*(row_i-1)+i} = row_nodes{row_i}{i};

                % Check if they're touching
                touch_check{i} = intersect(row_nodes{row_i}{i}, row_nodes{row_i-1}{i},"rows");
            end
        end
    end

    % Build the rows above this
    for row_i = middle_row - 1 : -1 : 1
        % Space the row out and find it
        max_vals         = max(cell2mat(cellfun(@(x) max(x), row_nodes{row_i + 1}, 'UniformOutput', false))); % find the bottom of the previous row
        [row_center, ~]  = Find_Center(local_nodes, max_vals(3) + dist, E, FoB);                              % find the closest center point after walking down the vert gap and half electrode height
        row_nodes{row_i} = create_rows(local_nodes, row_center, E, all_nodes, body_faces);                    % make the row
        
        % Set aside the electrode and check if rows are touching top-to-bottom
        touch_check = cell(1,E.E_count(2));
        for i = 1:E.E_count(2)
            % Save the individual electrode
            E_nodes{E.E_count(2)*(row_i-1)+i} = row_nodes{row_i}{i};

            % Check if they're touching
            touch_check{i} = intersect(row_nodes{row_i}{i}, row_nodes{row_i+1}{i},"rows");
        end

        % If touching, move further down
        search_dist = 0;
        while any(cellfun(@(x) ~isempty(x), touch_check))
            search_dist      = search_dist + 1; % Look 1 mm further away
            [row_center, ~]  = Find_Center(local_nodes, max_vals(3) + dist + search_dist, E, FoB); % find the closest center point after walking down the vert gap and half electrode height
            row_nodes{row_i} = create_rows(local_nodes, row_center, E, all_nodes, body_faces);     % make the row
            for i = 1:E.E_count(2)
                % Save the individual electrode
                E_nodes{E.E_count(2)*(row_i-1)+i} = row_nodes{row_i}{i};

                % Check if they're touching
                touch_check{i} = intersect(row_nodes{row_i}{i}, row_nodes{row_i+1}{i},"rows");
            end
        end
    end
end

function E_nodes = create_belts(local_nodes, E_plane, E, all_nodes, body_faces, flags, perim_mm)
    %{
    Create belts of electrodes and find the nodes/connectivity
    9/27/24 - Kyler Howard

    param: local_nodes - All nodes on to look at for making the electrode
    param: E_plane     - A subset of nodes to look through to find centers of electrodes quickly
    param: E           - Structure containing the size and type of electrode

    return: E_nodes   - Coordinates of the electrodes
    %}

    nheights = size(E.E_count,2);
    E_nodes  = cell(1,nheights);
    for i = 1:nheights
        E_nodes{i} = cell(E.E_count(i), 1);
    
        center = (max(E_plane{i},[],1) + min(E_plane{i},[],1)) / 2;
    
        if E.equal_space == 1 
    
            % KH: Equal Arc Length Electrodes
            ii = 1;
            jj = 1;
            arc_length = 0;
            n_points   = 2000;
            point      = zeros(n_points, 3);
        
            % Determine the order to place the electrodes
            if flags.use_GE == 1
                tht_i = 3*pi/2;
                d_tht = -(2*pi)/n_points;
                tht_f = -pi/2;
            elseif flags.use_GE == 0
                tht_i = pi;
                d_tht = 2*pi/n_points;
                tht_f = 3*pi;
            end
        
            for theta = tht_i : d_tht : tht_f - d_tht
                radius =  Parameratize_Bdry(E_plane{i}, 40, theta);
                point(jj,:) = [center(1) + radius*cos(theta), center(2) + radius*sin(theta), center(3)];
                
                % Make the first electrode
                if theta == tht_i
                    c_point       = point(jj,:);
                    E_center      = E_plane{i}(dsearchn(E_plane{i},c_point),:);
                    E_nodes{i}{ii} = create_electrode(local_nodes, E_center, E, all_nodes, body_faces);
                    ii = ii + 1;
                end
        
                % Measure the distance between sweeps
                if jj >= 2  
                    dist = norm(point(jj-1,1:2) - point(jj,1:2));
                    % dist = sqrt((point(j-1,1) - point(j,1))^2 + (point(j-1,2) - point(j,2))^2);
                    arc_length = arc_length + dist;
                end
        
                % Check if we have moved around enough
                % goal_arc_length = perim_mm{i} / (E.E_count(i) + 0.5);
                goal_arc_length = perim_mm{i} / (E.E_count(i) + 0);
                if arc_length >= goal_arc_length
                    % Reset arc length
                    arc_length = 0;
        
                    % Find the electrode
                    c_point       = point(jj,:);
                    E_center      = E_plane{i}(dsearchn(E_plane{i},c_point),:);
                    E_nodes{i}{ii} = create_electrode(local_nodes, E_center, E, all_nodes, body_faces);
                    ii = ii + 1;
        
                    % if i > E.E_count(i)
                    %     break
                    % end
                end
        
                % Update the point index
                jj = jj + 1;
            end
        
        else % Unequal spacing for GE
            % KH: One Quarter at a time effectively
            ii = 1;
            jj = 1;
            n_points   = 2000;
            point      = zeros(n_points, 3);
        
            % Determine the order to place the electrodes
            if flags.use_GE == 1
                tht_is = [pi, pi, 0, 0];
                d_tht = (2*pi)/n_points;
                tht_fs = [3*pi/2, pi/2, pi/2, -pi/2];
    
                thetas = [tht_is(1):d_tht:tht_fs(1)-d_tht, tht_is(2):-d_tht:tht_fs(2)+d_tht, tht_is(3):d_tht:tht_fs(3)-d_tht, tht_is(4):-d_tht:tht_fs(4)+d_tht];
            elseif flags.use_GE == 0
                error("Why are you using unequal spacing for ACT 5?")
            end
        
            for theta = thetas
    
                radius =  Parameratize_Bdry(E_plane{i}, 40, theta);
                point(jj,:) = [center(1) + radius*cos(theta), center(2) + radius*sin(theta), center(3)];
                
                if theta == tht_is(ceil(ii/4))
                    k = 1;
                    arc_length = 0;
                else % Measure the distance between sweeps 
                    dist = sqrt((point(jj-1,1) - point(jj,1))^2 + (point(jj-1,2) - point(jj,2))^2);
                    arc_length = arc_length + dist;
                end
        
                % Check if we have moved around enough. First electrode for each quadrant has to be shorter distance
                if mod(ii,4) == 1 % First electrode for each quadrant
                    if E.shape == "circle"
                        goal_arc_length = E.E_space / 2 + E.E_rad;
                    elseif E.shape == "rectangle"
                        goal_arc_length = E.E_space / 2 + E.E_width / 2;
                    end
                else % Remaning 3 electrodes
                    if E.shape == "circle"
                        goal_arc_length = E.E_space + E.E_dia;
                    elseif E.shape == "rectangle"
                        goal_arc_length = E.E_space + E.E_width;
                    end
                end
    
                if arc_length > goal_arc_length
                    % Reset arc length
                    arc_length = 0;
        
                    % Find the electrode
                    c_point    = point(jj,:);
                    E_center   = E_plane{i}(dsearchn(E_plane{i},c_point),:);
                    E_nodes{i}{ii} = create_electrode(local_nodes, E_center, E, all_nodes, body_faces);
                    ii = ii + 1;
                    k = k + 1;
        
                    if ii > E.E_count(i)
                        break
                    end
                end
    
                if k == 5 && mod(ii,4) == 1
                    arc_length = 0;
                end
        
                % Update the point index
                jj = jj + 1;
            end
    
            % Reorder the electrodes to be in the right order from the back, instead of going out from each side
            new_ind = [4:-1:1, 5:8, 12:-1:9, 13:16];
            E_nodes{i} = E_nodes{i}(new_ind);
        end
    end
end


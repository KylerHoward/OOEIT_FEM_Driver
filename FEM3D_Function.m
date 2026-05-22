function n_bframes = FEM3D_Function(filepath, filename, sbj_name, sbj_save_path, flags, noise)
    %{
    Run a 3D FEM simulation on the subject selected with the given settings
    The driver expects the subject to have the origin at the bottom, posterior,
    left corner of the torso. Positive y is towards the anterior side.
    All units NEED to be in standard units to work.
    Initial Version: 10/07/24 - Kyler Howard
    Current Version: 01/30/26 - Kyler Howard
    
    param: filepath - Filepath to the GMSH mesh
    param: filename - Filename of the GMSH mesh
    param: sbj_name - Character array of the subject's name
    param: flags    - Settings to run the simulation
    param: noise    - Noise and error vector

    load: tetmesh  - Mesh structure containing the nodes, connectivity, and labels
    load: cur_pat  - Current patterns for belt/patch electodes
    load: sbj_info - Structure containing heights and other info for each subject
    
    save: Umeas         - Measured voltages on each electrode for each current pattern (with or without noise appended to the name)
    save: Uall          - Measured voltages on all nodes for each current pattern (with or without noise appended to the name)
    save: cur_pat       - Current pattern used to create the simulation
    save: perim_mm      - Perimeter of the subject in mm
    save: noise         - Noise and error vector used to create the simulation
    save: sigma         - Conductivity on every node
    save: sigma_GT      - Ground truth conductivity matrix. One slice for each plane of electrodes
    save: nodes         - All nodes in space
    save: E_connect     - Connectivity matrix for each electrode
    save: flags         - Settings used to create the simulation
    save: volt_metadata - All metadata related to volt saving. Includes units and variable descriptions
    save: cond_metadata - All metadata related to cond saving. Includes units and variable descriptions
    %}
       
% ----------------------------------------------------------------------- %
%%                                 Setup                                  %
% ----------------------------------------------------------------------- %
    if flags.verbose == 1
        fprintf("   Defining Variables\n")
    end

    % Load the mesh & split the structure
    try % regular mesh from GMSH
        tetmesh      = load(fullfile(filepath,filename), "tetmesh").tetmesh;
        connectivity = double(tetmesh.cell.'); % indices of nodes making up an element
        nodes        = double(tetmesh.node.'); % xyz coordinates in mm
        labels       = tetmesh.field.';        % material id for each element
    catch % Saved heart meshes
        tetmesh      = load(fullfile(filepath,filename));
        connectivity = double(tetmesh.connectivity); % indices of nodes making up an element
        nodes        = double(tetmesh.nodes); % xyz coordinates in mm
        labels       = tetmesh.labels;        % material id for each element

        flags.are_bones = 0;
    end
    clear tetmesh
    
    % Setting the Injection Current Pattern
    if flags.inject_current == 1
        if flags.E_type == "patch"
            if flags.CP_choice == 1
                cur_pat = load(fullfile("Current_Patterns", "Box4by4aMxStrength8.mat"), "Box4by4aMxStrength8").Box4by4aMxStrength8;
            elseif flags.CP_choice == 2
                % current_pattern = load(fullfile("Current_Patterns", "CP32_4x8.mat"), "CP32_4x8").CP32_4x8;
                % Current pattern from GE is in uA
                cur_pat = load(fullfile("Current_Patterns", "GE_CP32_4x8.mat"), "CP32_4x8_Unscaled").CP32_4x8_Unscaled;
        
                % Scale the current pattern
                CP_scale        = load(fullfile("Current_Patterns", "GE_CP32_4x8.mat"), "CP32_4x8_IScale").CP32_4x8_IScale;
                CP_scale        = repmat(CP_scale.', 1, size(cur_pat,2));
                cur_pat = CP_scale .* cur_pat;
                
                % Convert to be in Amps
                cur_pat = cur_pat * 1e-6; 
            end
        elseif flags.E_type == "belt"
            % current_pat = load(fullfile("Current_Patterns", "CP32_16x2_M1.mat"), "Cur_pat3D").Cur_pat3D; % Normalized CP
            if flags.use_GE == 1
                % Current pattern from GE is in uA
                cur_pat = load(fullfile("Current_Patterns", "clinical_belt_CP.mat"), "CP_belt_16x2").CP_belt_16x2;
        
                % Scale the current pattern
                CP_scale        = load(fullfile("Current_Patterns", "clinical_belt_CP.mat"), "CPscale_belt_16x2").CPscale_belt_16x2;
                CP_scale        = repmat(CP_scale.', 1, size(cur_pat,2));
                cur_pat = CP_scale .* cur_pat;
                
                % Convert to be in Amps
                cur_pat = cur_pat * 1e-6;
            else
                % Current pattern from ACT5 is in Amps
                cur_pat = load(fullfile("Current_Patterns", "ACT5_CP32_2x16.mat"), "cur_pattern").cur_pattern;
                cur_pat = cur_pat(:,1:31);
            end
        end
    else % You only want to measure, not inject current
        if flags.E_choice <= 4
            L = 32;
        elseif flags.E_type == "patch"
            L = flags.E_count(1)*flags.E_count(2) * 2;
        elseif flags.E_type == "belt"
            L = flags.E_count * 2;
        end
        cur_pat = zeros(L, 1);
    end
    
    L = size(cur_pat, 1); % Number of electrodes
    K = size(cur_pat, 2); % Number of current patterns
    % if K == 1 % KH 1/28/26: I'm not sure why it was set to 0
    %     K = 0;
    % end
    
    % Load the table
    if contains(sbj_name, "R1")
        sbj_sheet = readtable("CT Data Boundaries.xlsx", "Sheet","Lungmap Set");
    elseif contains(sbj_name, "case1")
        sbj_sheet = readtable("CT Data Boundaries.xlsx", "Sheet","New Mexico Baby Set");
    elseif contains(sbj_name, "EIT1")
        sbj_sheet = readtable("CT Data Boundaries.xlsx", "Sheet","Anschutz Set");
    elseif contains(sbj_name, "EIT2")
        sbj_sheet = readtable("CT Data Boundaries.xlsx", "Sheet","R01");
    elseif contains(sbj_name, "G2")
        sbj_sheet = readtable("CT Data Boundaries.xlsx", "Sheet","R21");
    elseif contains(sbj_name, "MCR0")
        sbj_sheet = readtable("CT Data Boundaries.xlsx", "Sheet","MCR Set");
    else
        error("No excel file for subject %s\n", sbj_name)
    end
    
    % Extract wanted info from the table
    row_num       = find(strcmp(sbj_sheet.Subject, sbj_name));
    carina_height = sbj_sheet.CarinaFromST_mm(row_num);
    T5_height     = sbj_sheet.T5FromST_mm(row_num);
    T8_height     = sbj_sheet.T8FromST_mm(row_num);
    
    if flags.are_bones == 1
        background  = 1;
        lung        = 2;
        trachea     = 3;
        soft_tissue = 4;
        bone        = 5;
        esophagus   = 6;
        heart       = 7;
        external    = 8;
    else
        lung        = 1;
        trachea     = 2;
        soft_tissue = 3;
        esophagus   = 4;
        heart       = 5;
        background  = 6;
        external    = 7;
    end
    
    % Create a cell array with each individual organ mesh
    organ_connects = cell(length(unique(labels)),1);
    i = 1;
    for label = unique(labels).'
        organ_connects{i} = connectivity(labels==label, :);
        i = i + 1;
    end
    clear i
    
% ----------------------------------------------------------------------- %
%%                        Rotation and Translation                        %
% ----------------------------------------------------------------------- %
    % Rotate and then translate the body to standard axis:
        % Anterior is positive y
        % Superior is positive z
        % Right    is positive x
    [nodes, sbj_info] = Rotate_and_Translate_Body(nodes, organ_connects, carina_height, T5_height, T8_height, flags);
    
% ----------------------------------------------------------------------- %
%%                              Surface Nodes                             %
% ----------------------------------------------------------------------- %
    if flags.verbose == 1
        fprintf("   Extracting Trachea and Surface Nodes\n")
    end

    % Extract desired faces & find their intersection
    body_faces      = Get_Unique_Faces(organ_connects{soft_tissue});
    lung_faces      = Get_Unique_Faces(organ_connects{lung});
    trachea_faces   = Get_Unique_Faces(organ_connects{trachea});
    heart_faces     = Get_Unique_Faces(organ_connects{heart});
    esophagus_faces = Get_Unique_Faces(organ_connects{esophagus});
    
    lung_intersect      = intersect(sort(body_faces,2),  sort(lung_faces,2),      "rows");
    trachea_intersect   = intersect(sort(body_faces,2),  sort(trachea_faces,2),   "rows");
    heart_intersect     = intersect(sort(body_faces,2),  sort(heart_faces,2),     "rows");
    esophagus_intersect = intersect(sort(body_faces,2),  sort(esophagus_faces,2), "rows");
     
    inner_intersects = vertcat(lung_intersect, trachea_intersect, heart_intersect, esophagus_intersect);
    surface_faces    = setdiff(sort(body_faces,2), sort(inner_intersects,2), "rows");
    
    % Extract desired organ nodes
    [trachea_nodes, ~]    = Get_Tet_Nodes(nodes, organ_connects{trachea});
    [lung_nodes, ~]       = Get_Tet_Nodes(nodes, organ_connects{lung});
    [boundary_nodes, ind] = Get_Surface_Nodes(nodes, surface_faces);
    
    % Delete any extra nodes in the center of boundary_nodes
    inside_indices                   = Find_Internal_Nodes(boundary_nodes, flags);
    old_boundary_nodes               = boundary_nodes;
    boundary_nodes(inside_indices,:) = [];
    ind(inside_indices,:)            = [];
    
    if flags.plot_trachea == 1
        figure()
        subplot(1,2,1)
            scatter3(trachea_nodes(:,1), trachea_nodes(:,2), trachea_nodes(:,3), "MarkerEdgeAlpha", 0.2)
            hold on
            scatter3(mean(trachea_nodes(:,1)), mean(trachea_nodes(:,2)), sbj_info.carina, "r", "filled")
            scatter3(trachea_nodes(startsWith(string(trachea_nodes(:,3)), sprintf("%.1f", sbj_info.carina)), 1), trachea_nodes(startsWith(string(trachea_nodes(:,3)), sprintf("%.1f", sbj_info.carina)), 2), sbj_info.carina, "r", "filled")
            % constantplane("z", sbj_info.carina) R2024b
            title(sprintf("Carina: %.2f mm", sbj_info.carina))
            axis equal
        subplot(1,2,2)
            pdeplot3D(nodes.', organ_connects{lung}.')
    end
    
    if flags.save_heart_mesh == 1
        [heart_nodes, ~]         = Get_Tet_Nodes(nodes, organ_connects{heart});
        [heart_surface_nodes, ~] = Get_Surface_Nodes(nodes, heart_faces);
    
        heart_name = sprintf("%s_Heart_Mesh.mat", sbj_name);
        save(fullfile("Heart_Meshes", heart_name), "nodes", "connectivity", "labels", "heart_faces", "heart_nodes", "heart_surface_nodes")
    end
    
    clear trachea_faces heart_faces esophagus_faces
    clear lung_intersect trachea_intersect heart_intersect esophagus_intersect inner_intersects surface_faces
% ----------------------------------------------------------------------- %
%%                             Make Electrodes                            %
% ----------------------------------------------------------------------- %
    if flags.verbose == 1
        fprintf("   Making Electrodes\n")
    end

    tic
    [E_nodes, perim_mm] = Make_Electrodes3(boundary_nodes, nodes, body_faces, sbj_info, flags);
    if iscell(E_nodes) == 0
        E_nodes = {E_nodes};
    end
    
    % Trim electrodes that touch
    [E_nodes, percent_trimmed] = Trim_Electrodes(E_nodes, flags);

    % Electrodes are too big and we have a smaller belt/patch available
    if sum(flags.E_choice == [1,3]) && percent_trimmed > 2
        if flags.verbose == 1
            fprintf("      Electrodes are too big. Trying again at smaller size\n")
        end
        flags.E_choice = flags.E_choice + 1;

        [E_nodes, perim_mm] = Make_Electrodes3(boundary_nodes, nodes, body_faces, sbj_info, flags);
        if iscell(E_nodes) == 0
            E_nodes = {E_nodes};
        end
        
        % Trim electrodes that touch
        [E_nodes, ~] = Trim_Electrodes(E_nodes, flags);
    end

    % Find the surface mesh faces that make up the electrodes
    E_connect = Align_Electrode_Faces(nodes, body_faces, E_nodes, flags);
    
    make_time = toc;
    if flags.verbose == 1
        fprintf("      It took %.2f seconds to make the electrodes\n", make_time)
        fprintf("      The perimeter is %.2f mm\n", perim_mm)
    end
    
    % Look at the electrodes and plot their number of faces & areas per electrode
    Compare_Electrodes(L, nodes, E_connect, flags)
    
    
    if flags.plot_electrodes == 1
        figure()
            hold on
            % scatter3(boundary_nodes(:,1), boundary_nodes(:,2), boundary_nodes(:,3), "MarkerEdgeColor", [0.3010 0.7450 0.9330], "LineWidth",0.05)
            scatter3(boundary_nodes(:,1), boundary_nodes(:,2), boundary_nodes(:,3), "MarkerEdgeColor", [0.21 0.71 0.52], "LineWidth",0.05)
            scatter3(trachea_nodes(:,1),trachea_nodes(:,2),trachea_nodes(:,3),"y","filled")
            scatter3(lung_nodes(:,1), lung_nodes(:,2), lung_nodes(:,3), "b")
            for i = 1:length(E_connect)
                trimesh(E_connect{i}, nodes(:,1), nodes(:,2), nodes(:,3), "FaceColor", "r", "edgeColor", "r");
            end
            xlabel("X (mm)")
            ylabel("Y (mm)")
            zlabel("Z (mm)")
            axis equal
        if flags.do_pauses == 1
            fprintf("      Press any key if correct\n")
            pause
        end
    end
    
% ----------------------------------------------------------------------- %
%%                            Heart Conditions                            %
% ----------------------------------------------------------------------- %
    % Setting up heart boundary condition information
    if flags.heart_BCs == 1
        if flags.verbose == 1
            fprintf("   Creating heart Dirchlet boundary conditions\n")
            heart_start_time = tic;
        end
        Diriclet_name = dir(fullfile("Heart_BCs", sprintf("%s_*Dirichlet*.mat",sbj_name))).name;
        heart_BCs = load(fullfile("Heart_BCs", Diriclet_name));

        % Find Heart unit conversion
        try
            % Find the metadata substructure and go straigt to the units
            heart_var_names = fieldnames(heart_BCs);
            meta_idx = contains(lower(heart_var_names), "metadata");
            heart_metadata = heart_BCs.(heart_var_names{meta_idx}).units;
            
            % Find the unit used
            meta_var_names = fieldnames(heart_metadata);
            unit_idx = contains(lower(meta_var_names), "values");
            unit_name = heart_metadata.(meta_var_names{unit_idx});
            if contains(unit_name, "m")
                heart_unit_convert = 1e-3;
            elseif contains(unit_name, "u") || contains(unit_name, "µ")
                heart_unit_convert = 1e-6;
            else
                heart_unit_convert = 1;
            end
        catch
            heart_unit_convert = 1e-3;
        end
        
        heart_BC_vals            = heart_BCs.heart_BC_values*heart_unit_convert; % Convert to standard units
        heart_BC_surface_indices = heart_BCs.heart_BC_surface_indices;
        clear heart_BCs
    
        % Convert the given indices into the node coordinates
        heart_surface_nodes   = load(fullfile("Heart_Meshes", sprintf("%s_Heart_Mesh.mat",sbj_name))).heart_surface_nodes;
        % heart_faces           = load(fullfile("Heart_Meshes", sprintf("%s_Heart_Mesh.mat",sbj_name))).heart_faces;

        % Find the heart nodes that we're applying the BC on
        n_hframes        = size(heart_BC_surface_indices, 1);
        heart_BC_indices = cell(size(heart_BC_surface_indices));
            
        for hframe = 1:n_hframes
            for i_BC = 1:size(heart_BC_surface_indices,2)
                heart_BC_nodes = heart_surface_nodes(heart_BC_surface_indices{hframe,i_BC},:);
                
                % Find the global nodes indicies for the heart BC
                [~, heart_BC_indices{hframe,i_BC}] = intersect(nodes, heart_BC_nodes, "rows", "stable");
            end
        end
        if flags.plot_heart == 1
            % Set figure settings
            fig = uifigure("color","w","Position",[573,337.67,560,470]);
            % --- Grid Layout ---
            % Layout:
            %   Row 1: plots (1 or 2)
            %   Row 2: sliders (3 sliders stacked) autofitted
            mainGrid = uigridlayout(fig,[2 1]);
            mainGrid.RowHeight = {"1x", "fit"};
            mainGrid.ColumnWidth = {"1x"};
    
            plotGrid = uigridlayout(mainGrid,[1 1]);
            ax1 = uiaxes(plotGrid);
            hold(ax1, "on")
            start_hframe = 5;
            
            % Plot 1 for Frame 1
            s1 = cell(1,size(heart_BC_surface_indices,2));
            for i_BC = 1:size(heart_BC_surface_indices,2)
                heart_BC_nodes = heart_surface_nodes(heart_BC_surface_indices{start_hframe,i_BC},:);
                s1{i_BC} = scatter3(ax1, heart_BC_nodes(:,1), heart_BC_nodes(:,2), heart_BC_nodes(:,3),[],heart_BC_vals(i_BC)*ones(size(heart_BC_nodes,1),1),"filled");
            end
            title(ax1, sprintf("Heart BC\nFrame %d", 1))
            xlabel(ax1, "X (mm)")
            ylabel(ax1, "Y (mm)")
            zlabel(ax1, "Z (mm)")
            axis(ax1, "equal")
            clim(ax1, [min(heart_BC_vals), max(heart_BC_vals)])
            colormap(ax1, "turbo")
            c = colorbar(ax1, "eastoutside");
            c.Label.String = "Volts";

            % --- Frame Slider Block --- 
            frameBlock = uigridlayout(mainGrid,[2 1]); 
            frameBlock.RowHeight = {15,30};
            
            uilabel(frameBlock, ... 
                    "Text","Frame Index", ... 
                    "HorizontalAlignment","center", ... 
                    "FontWeight","bold");
            uislider(frameBlock, ...
                     "Limits",[1 n_hframes], ...
                     "Value",start_hframe, ...
                     "MajorTicks",1:n_hframes, ...
                     "MinorTicks",[], ...
                     "ValueChangedFcn",@(src,evt) updateHeartFrame(ax1, src.Value, heart_BC_vals, heart_surface_nodes, heart_BC_surface_indices));
        end

        if flags.verbose == 1
            heart_stop_time = toc(heart_start_time);
            fprintf("      It took %.2f seconds to set heart Dirichlet conditions\n", heart_stop_time)
        end
    else
        % Set empty conditions for no Dirichlet BCs
        heart_BC_indices = {[]};
        heart_BC_vals  = [];
        n_hframes = 1;
    end
    
% ----------------------------------------------------------------------- %
%%                             Create Solver                              %
% ----------------------------------------------------------------------- %
    if flags.solve_problem == 1
        % Create a mesh-object from nodes, connections, and electrode connections (NODES IN METERS)
        if flags.verbose == 1
            fprintf("   Making Forward Mesh\n")
        end
        if flags.do_parfor == 1
            fmesh = ParForwardMesh1st(nodes*1e-3, connectivity, E_connect);
        else
            fmesh = ForwardMesh1st(nodes*1e-3, connectivity, E_connect);
        end
        
        % Initialize the forward problem solver (FPS). This is the object, that
        % (utilizing the mesh-object given to it) computes forward problem solutions
        % using the FEM approximation of the complete electrode model (CEM).
        if flags.verbose == 1
            fprintf("   Initializing FEM Solver\n")
        end
        init_start = tic;
        if flags.do_parfor == 1
            solver = MF_EITFEM(fmesh);
        else
            solver = EITFEM(fmesh);
        end
        init_time  = toc(init_start);
        if flags.verbose == 1
            fprintf("      It took %.2f minutes to initialize solver\n", init_time/60)
        end
                
        % Set solving mode to injecting current and measuring voltages
        % Set the contact impedance - had 2.4 originally in 2D FEM, then 4.8, then 0.05
        % Fine testing with closly matched rect. electrodes:  0.2475 for Subj011
        % Fine testing with closly matched patch. electrodes: 0.05 for Kyler
        % flags.zeta = 0.05*ones(L,1); % R1133
        % flags.zeta = 0.116*ones(L,1); % R1044 for sbj005
        % flags.zeta = 0.079*ones(L,1); % case126002 for sbj002
    
        % Creating an equation based on data for sbj002 & sbj005. (339.62, 0.079) & (290.67, 0.116)
        % FIXME: KH 1/26/26, this is only for babies on GE right now
        % zeta       = round((0.116-0.079)/(290.67-339.62) * (perim_mm - 290.67) + 0.116,3);

        % Creating zeta based on the average of the first four GE subjects
        % FIXME: KH 5/21/26, this is only for babies on GE right now
        zeta = mean([0.079, 0.116, 0.217, 0.133]);
        flags.zeta = zeta*ones(L,1);
        solver.mode = "current";
        solver.zeta = flags.zeta;
    end

    % % DELETE ME: TESTING CONTACT IMPEDANCES
    % flags2 = flags;
    % solver2 = solver;
    % zeta = [0.082 0.083 0.084 0.086 0.087 0.088 0.089];
    % parfor ii = 1:length(zeta)
    % 
    %     % Create local copies of the flags and the solver
    %     flags  = flags2;
    %     solver = solver2;
    % 
    %     solver.zeta = zeta(ii)*ones(L,1);
    %     flags.zeta  = solver.zeta;  
    
% ----------------------------------------------------------------------- %
%%                               Setup Loops                              %
% ----------------------------------------------------------------------- %
    if flags.do_conditions == 1
        num_conditions = numel(flags.conditions);
    else
        num_conditions = 1;
        condition_name = "Custom";
    end
    for i_condition = 1:num_conditions
        if flags.do_conditions == 1
            % Save the condition settings
            condition = flags.conditions{i_condition};
            flags.max_inspiration = condition{1};
            flags.equal_vent      = condition{2};
            flags.left_only       = condition{3};
            flags.right_only      = condition{4};
            flags.esoph_intubate  = condition{5};
            condition_name        = condition{6};
        end
            
        if flags.permute_conds == 1 && flags.do_conditions == 1
            % Save the permutation settings but don't overwrite if we are
            % setting just one permutation with the custom "condition"
            permutation       = flags.permutations{i_condition};
            num_permutations  = permutation{1};
            flags.lung_range  = permutation{2};
            flags.esoph_range = permutation{3};
        else
            num_permutations = 1;
        end

        for i_permutation = 1:num_permutations
            
            if flags.make_video == 1
                % Find how many frames we are simulating
                breaths_per_sec = 1 / (flags.breath_rate/60);
                n_bframes = ceil(breaths_per_sec * flags.fps);

                % Create cosine curve for ventilation
                cosA = range(flags.insp_range)/2;
                cosB = 2*pi / n_bframes;
                cosC = n_bframes / 2; 
                cosD = mean(flags.insp_range);

                flags.breath_curve = cosA * cos(cosB *((1:n_bframes) - cosC)) + cosD; 

                % Create cardiac cycle based on literature
                beats_per_sec = 1 / (flags.heart_rate/60);
                heart_period = ceil(beats_per_sec * flags.fps);
                flags.heart_curve = Make_Cardiac_Curve(heart_period, n_bframes);

                if flags.plot_breath == 1
                    figure
                        hold on
                        plot(1:n_bframes, flags.breath_curve)
                        plot(1:n_bframes, flags.heart_curve)
                        xline([2,20,37],"k")
                        legend(["Breath Curve", "Heart Curve","","",""],"Location","southoutside")
                        xlabel("Frame")
                        ylabel("% of Cycle")
                end
            else
                % Set the default values
                n_bframes = 1;
                flags.breath_curve = flags.max_inspiration;
                flags.heart_curve  = flags.cardiac_cycle;
            end

% ----------------------------------------------------------------------- %
%%                         Assign Conductivities                          %
% ----------------------------------------------------------------------- %
            % Initialize Umeas and Uall
            sigma = zeros(size(nodes,1), n_bframes);
            Umeas = zeros(L,             K, n_hframes*n_bframes);
            Uall  = zeros(size(nodes,1), K, n_hframes*n_bframes);
            Umeas_cell = cell(n_bframes, 1);
            Uall_cell  = cell(n_bframes, 1);

            % KH: Testing making this a parfor. Probably will need to check for computer availability?
            parfor bframe = 1:n_bframes
            % for bframe = 1:n_bframes
                % Create local copies of the flags and the solver
                flags_local  = flags;
                solver_local = solver;

                % Extract which point on the breath curve we want to simulate on a local copy of the flags
                flags_local.max_inspiration = flags.breath_curve(bframe);
                flags_local.cardiac_cycle   = flags.heart_curve(bframe);

                if flags_local.verbose == 1
                    fprintf("   Assigning Conductivities\n")
                end
            
                % Creating the conductivity vector at the nodes (IN SIEMENS PER METER)
                sigma(:,bframe) = Assign_Conductivities(nodes, connectivity, labels, lung_nodes, flags_local);
            
    % ----------------------------------------------------------------------- %
    %%                                 Solve                                  %
    % ----------------------------------------------------------------------- %                
                % Make a local version for each core
                Umeas_local = zeros(L, K, n_hframes);
                Uall_local  = zeros(size(nodes,1), K, n_hframes);

                for hframe = 1:n_hframes
                    if flags_local.solve_problem == 1
                        fprintf("   Solving Forward Problem %d of %d\n", (i_permutation-1)*n_hframes*n_bframes + (bframe-1)*n_hframes + hframe, num_permutations*n_hframes*n_bframes)
                    
                        % Set the current pattern (IN AMPS)
                        solver_local.Iel = cur_pat;
                    
                        % Solve the voltage
                        solve_start          = tic;
                        if flags_local.do_parfor == 1
                            [Umeas_frame, Imeas_frame, Uall_frame] = MF_Simulation(fmesh,  [], sigma(:,bframe), solver_local.zeta, heart_BC_indices(hframe,:), heart_BC_vals, solver, "current", noise);
                        else
                            [Umeas_frame, Imeas_frame, Uall_frame] = MF_Simulation2(fmesh, [], sigma(:,bframe), solver_local.zeta, heart_BC_indices(hframe,:), heart_BC_vals, solver, "current", noise);
                        end
                        solve_time           = toc(solve_start);
                        if flags_local.verbose == 1
                            fprintf("      It took %.2f minutes to solve for voltages\n", solve_time/60)
                        end
                    
                        Umeas_frame = Umeas_frame * 1e3; % Convert V to mV
                        Uall_frame  = Uall_frame  * 1e3; % Convert V to mV

                        % Save the local frame
                        Umeas_local(:,:,hframe) = reshape(Umeas_frame, L, K); 
                        Uall_local(:,:,hframe) = Uall_frame(1:size(nodes,1),:);
                        % Umeas(:,:,(bframe-1)*n_hframes + hframe) = reshape(Umeas_frame, L, K);
                        % Uall(:,:,(bframe-1)*n_hframes + hframe)  = Uall_frame(1:size(nodes,1),:);
                    end % end running a solve
                end % end looping over heart frames

                % Combine all the heart solutions into a cell array
                Umeas_cell{bframe} = Umeas_local; 
                Uall_cell{bframe} = Uall_local;

            end % end looping over each video

            % Finally combine all results
            Umeas = cat(3, Umeas_cell{:});
            Uall  = cat(3, Uall_cell{:});
                    
    % ----------------------------------------------------------------------- %
    %%                                Plotting                                %
    % ----------------------------------------------------------------------- %
            % Plot global nodal voltages
            if flags.solve_problem == 1 && flags.plot_volts > 0
                Plot_Voltages(nodes, Uall, flags)   
            
                if flags.do_pauses == 1
                    fprintf("      Press any key if correct\n")
                    pause
                end
            end
            
            % Plot the ground truth images at each row
            GT_Thickness = 5; %mm
            sigma_GT = Create_GT_Images(GT_Thickness, nodes, E_nodes, sigma, flags);
% ----------------------------------------------------------------------- %
%%                               Save Data                                %
% ----------------------------------------------------------------------- %
            
            if flags.solve_problem == 1

                if noise(1) == 0 || noise(2) == 0
                    Umeas_NoNoise = Umeas;
                    Uall_NoNoise  = Uall;
                end
            
                % Create the suffix for saving the file based on the settings chosen
                save_suffix = Make_Save_Name(condition_name, i_permutation, solver.zeta, flags);
        
                % Create metadata to include with the saving
                volt_metadata = Make_Metadata("Volt");
                cond_metadata = Make_Metadata("Cond");
            
                % Check if the user isn't just saving to the results folder
                if contains(sbj_save_path, fullfile("OOEIT_FEM_Driver","Results")) ~= 1
                    cond_save_path = fullfile(sbj_save_path, flags.E_type, "conductivities", condition_name);
                    volt_save_path = fullfile(sbj_save_path, flags.E_type, "voltages",       condition_name);
                else
                    cond_save_path = sbj_save_path;
                    volt_save_path = sbj_save_path;
                end
    
                % Save the voltages and conductivties
                volt_name = sprintf("%s-Volt%s.mat", sbj_name, save_suffix);
                cond_name = sprintf("%s-Cond%s.mat", sbj_name, save_suffix);
                % if noise(1) == 0 || noise(2) == 0
                %     save(fullfile(volt_save_path, volt_name), "Umeas_NoNoise", "Uall_NoNoise", "cur_pat", "perim_mm", "noise", "flags", "volt_metadata", "-v7.3")
                % else
                %     save(fullfile(volt_save_path,volt_name), "Umeas", "Uall", "cur_pat", "perim_mm", "noise", "flags", "volt_metadata", "-v7.3")
                % end
                % save(fullfile(cond_save_path,cond_name), "sigma", "sigma_GT", "nodes", "E_connect", "flags", "cond_metadata", "-v7.3")
                saveData(volt_save_path, cond_save_path, volt_name, cond_name,  Umeas_NoNoise, Uall_NoNoise, Umeas, Uall, cur_pat, perim_mm, noise, flags, volt_metadata, sigma, sigma_GT, nodes, E_connect, cond_metadata)
            end % end saving data
        end % end looping over permutations
    end % end looping through conditions
    % end % Zeta testing end
end % end function as a whole
    
% ----------------------------------------------------------------------- %
%%                            Custom Functions                            %
% ----------------------------------------------------------------------- %

function updateHeartFrame(ax1, framenum, heart_BC_vals, heart_surface_nodes, heart_BC_surface_indices)
    plotframe = round(framenum);
    cla(ax1)
    hold(ax1,"on")

    % Update plot 1
    s1 = cell(1,size(heart_BC_surface_indices,2));
    for i_BC = 1:size(heart_BC_surface_indices,2)
        heart_BC_nodes = heart_surface_nodes(heart_BC_surface_indices{plotframe,i_BC},:);
        colordata = heart_BC_vals(i_BC)*ones(size(heart_BC_nodes,1),1);
        s1{i_BC} = scatter3(ax1, heart_BC_nodes(:,1), heart_BC_nodes(:,2), heart_BC_nodes(:,3),[],colordata,"filled");
    end
    title(ax1, sprintf("Heart BC\nFrame %d", 1))
    xlabel(ax1, "X (mm)")
    ylabel(ax1, "Y (mm)")
    zlabel(ax1, "Z (mm)")
    axis(ax1, "equal")
    clim(ax1, [min(heart_BC_vals), max(heart_BC_vals)])
    colormap(ax1, "turbo")
    c = colorbar(ax1, "eastoutside");
    c.Label.String = "Volts";

    
    title(ax1, sprintf("Heart BC\nFrame %d", plotframe))
end


function saveData(volt_save_path, cond_save_path, volt_name, cond_name, Umeas_NoNoise, Uall_NoNoise, Umeas, Uall, cur_pat, perim_mm, noise, flags, volt_metadata, sigma, sigma_GT, nodes, E_connect, cond_metadata)
    if noise(1) == 0 || noise(2) == 0
        save(fullfile(volt_save_path, volt_name), "Umeas_NoNoise", "Uall_NoNoise", "cur_pat", "perim_mm", "noise", "flags", "volt_metadata", "-v7.3")
    else
        save(fullfile(volt_save_path,volt_name), "Umeas", "Uall", "cur_pat", "perim_mm", "noise", "flags", "volt_metadata", "-v7.3")
    end
    save(fullfile(cond_save_path,cond_name), "sigma", "sigma_GT", "nodes", "E_connect", "flags", "cond_metadata", "-v7.3")
end
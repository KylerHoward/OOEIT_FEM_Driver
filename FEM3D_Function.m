function FEM3D_Function(filepath, filename, sbj_name, sbj_save_path, flags, noise)
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
    tetmesh      = load(fullfile(filepath,filename), "tetmesh").tetmesh;
    connectivity = double(tetmesh.cell'); % indices of nodes making up an element
    nodes        = double(tetmesh.node'); % xyz coordinates in mm
    labels       = tetmesh.field';        % material id for each element
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
                CP_scale        = repmat(CP_scale', 1, size(cur_pat,2));
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
                CP_scale        = repmat(CP_scale', 1, size(cur_pat,2));
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
    for label = unique(labels)'
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
    
    lung_intersect      = intersect(sort(body_faces,2),  sort(lung_faces,2),      'rows');
    trachea_intersect   = intersect(sort(body_faces,2),  sort(trachea_faces,2),   'rows');
    heart_intersect     = intersect(sort(body_faces,2),  sort(heart_faces,2),     'rows');
    esophagus_intersect = intersect(sort(body_faces,2),  sort(esophagus_faces,2), 'rows');
     
    inner_intersects = vertcat(lung_intersect, trachea_intersect, heart_intersect, esophagus_intersect);
    surface_faces    = setdiff(sort(body_faces,2), sort(inner_intersects,2), 'rows');
    
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
            scatter3(mean(trachea_nodes(:,1)), mean(trachea_nodes(:,2)), sbj_info.carina, 'r', 'filled')
            scatter3(trachea_nodes(startsWith(string(trachea_nodes(:,3)), sprintf("%.1f", sbj_info.carina)), 1), trachea_nodes(startsWith(string(trachea_nodes(:,3)), sprintf("%.1f", sbj_info.carina)), 2), sbj_info.carina, 'r', 'filled')
            % constantplane('z', sbj_info.carina) R2024b
            title(sprintf("Carina: %.2f mm", sbj_info.carina))
            axis equal
        subplot(1,2,2)
            pdeplot3D(nodes', organ_connects{lung}')
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
    E_nodes = Trim_Electrodes(E_nodes);
    
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
            scatter3(boundary_nodes(:,1), boundary_nodes(:,2), boundary_nodes(:,3), 'MarkerEdgeColor', [0.3010 0.7450 0.9330], 'LineWidth',0.05)
            scatter3(trachea_nodes(:,1),trachea_nodes(:,2),trachea_nodes(:,3),'y','filled')
            scatter3(lung_nodes(:,1), lung_nodes(:,2), lung_nodes(:,3), 'b')
            for i = 1:length(E_connect)
                trimesh(E_connect{i}, nodes(:,1), nodes(:,2), nodes(:,3), 'FaceColor', 'r', 'edgeColor', 'r');
            end
            xlabel('X (mm)')
            ylabel('Y (mm)')
            zlabel('Z (mm)')
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
        end
        heart_BCs = load(fullfile("Heart_BCs", sprintf("%s_Heart_BC_frame_5.mat",sbj_name)));
        
        heart_BC_currents     = heart_BCs.current_density*10^-3; % CONVERT mAMPS TO AMPS
        heart_BC_face_indices = heart_BCs.nodeGroups;
        clear heart_BCs

        heart_BCs = load(fullfile("Heart_BCs", sprintf("%s_Heart_BC_frame_10.mat",sbj_name)));
        heart_BC_face_indices = [heart_BC_face_indices; heart_BCs.nodeGroups];
    
        % Convert the given indices into the node coordinates
        heart_surface_nodes   = load(fullfile("Heart_Meshes", sprintf("%s_Heart_Mesh.mat",sbj_name))).heart_surface_nodes;
        heart_faces           = load(fullfile("Heart_Meshes", sprintf("%s_Heart_Mesh.mat",sbj_name))).heart_faces;

        % Find the heart nodes that we're applying the BC on
        n_hframes = size(heart_BC_face_indices, 1);
        heart_BC_node_indices = cell(size(heart_BC_face_indices));
        for hframe = 1:n_hframes
            if flags.plot_heart == 1
                figure;
            end
            for i_BC = 1:size(heart_BC_face_indices,2)
                heart_BC_faces                     = heart_faces(heart_BC_face_indices{hframe,i_BC},:);
                heart_BC_node_indices{hframe,i_BC} = unique(heart_BC_faces(:));
                heart_BC_nodes                     = nodes(heart_BC_node_indices{hframe,i_BC},:);
    
                if flags.plot_heart == 1
                    subplot(1,size(heart_BC_face_indices,2),i_BC)
                        hold on
                        scatter3(heart_surface_nodes(:,1), heart_surface_nodes(:,2), heart_surface_nodes(:,3),'k')
                        scatter3(heart_BC_nodes(:,1),      heart_BC_nodes(:,2),      heart_BC_nodes(:,3),'r','filled')
                        if i_BC == round(median(1:size(heart_BC_face_indices,2)))
                            title(sprintf("Heart Frame %d\n%.4f BC Condition", hframe, heart_BC_currents(i_BC)))
                        else
                            title(sprintf("\n%.4f BC Condition", heart_BC_currents(i_BC)))
                        end
                end
            end
        end
    else
        % Set empty conditions for no Dirichlet BCs
        heart_BC_node_indices = {[]};
        heart_BC_currents     = [];
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
        zeta       = round((0.116-0.079)/(290.67-339.62) * (perim_mm - 290.67) + 0.116,3);
        flags.zeta = zeta*ones(L,1);
        solver.mode = 'current';
        solver.zeta = flags.zeta;
    end
    
% ----------------------------------------------------------------------- %
%%                         Assign Conductivities                          %
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
            
        if flags.permute_conds == 1
            % Save the permutation settings
            permutation       = flags.permutations{i_condition};
            num_permutations  = permutation{1};
            flags.lung_range  = permutation{2};
            flags.esoph_range = permutation{3};
        else
            num_permutations = 1;
        end

        for i_permutation = 1:num_permutations
            
            if flags.verbose == 1
                fprintf("   Assigning Conductivities\n")
            end
        
            % Creating the conductivity vector at the nodes (IN SIEMENS PER METER)
            sigma = Assign_Conductivities(nodes, connectivity, labels, lung_nodes, flags);
        
    % ----------------------------------------------------------------------- %
    %%                       Create Ground Truth Images                       %
    % ----------------------------------------------------------------------- %
            % Plot the ground truth images at each row
            GT_Thickness = 5; %mm
            sigma_GT = Create_GT_Images(GT_Thickness, nodes, E_nodes, sigma, flags);
        
    % ----------------------------------------------------------------------- %
    %%                                 Solve                                  %
    % ----------------------------------------------------------------------- %
        % DELETE ME: TESTING CONTACT IMPEDANCES
        % zeta = [0.041];
        % for ii = 1:length(zeta)
        %     solver.zeta = zeta(ii)*ones(L,1);
        %     flags.zeta  = solver.zeta;    
            for hframe = 1:n_hframes
                if flags.solve_problem == 1
                    fprintf("   Solving Forward Problem %d of %d\n", (i_permutation-1)*n_hframes + hframe, num_permutations*n_hframes)
                
                    % Set the current pattern (IN AMPS)
                    solver.Iel = cur_pat;
                
                    % Solve the voltage
                    solve_start          = tic;
                    if flags.do_parfor == 1
                        [Umeas, Imeas, Uall] = MF_Simulation(fmesh, [], sigma, solver.zeta, heart_BC_node_indices, heart_BC_currents, solver, 'current', noise);
                    else
                        [Umeas, Imeas, Uall] = MF_Simulation2(fmesh, [], sigma, solver.zeta, heart_BC_node_indices, heart_BC_currents, solver, 'current', noise);
                    end
                    solve_time           = toc(solve_start);
                    if flags.verbose == 1
                        fprintf("      It took %.2f minutes to solve for voltages\n", solve_time/60)
                    end
                
                    Umeas = Umeas * 1e3; % Convert V to mV
                    Uall  = Uall  * 1e3; % Convert V to mV
                    Umeas = reshape(Umeas, L, K);
            
                    if noise(1) == 0 || noise(2) == 0
                        Umeas_NoNoise = Umeas;
                        Uall_NoNoise  = Uall;
                    end
                
                    % Create the suffix for saving the file based on the settings chosen
                    save_suffix = Make_Save_Name(condition_name, i_permutation, solver.zeta, hframe, flags);
            
                    % Create metadata to include with the saving
                    volt_metadata = Make_Metadata("Volt");
                    cond_metadata = Make_Metadata("Cond");
                
                    % Check if the user isn't just saving to the results folder
                    if contains(sbj_save_path, fullfile("OOEIT_FEM_Driver","Results")) ~= 1
                        cond_save_path = fullfile(sbj_save_path, flags.E_type, "conductivities", condition_name);
                        volt_save_path = fullfile(sbj_save_path, flags.E_type, "voltages",        condition_name);
                    else
                        cond_save_path = sbj_save_path;
                        volt_save_path = sbj_save_path;
                    end
    
                    % Save the voltages and conductivties
                    volt_name = sprintf("%s-Volt%s.mat", sbj_name, save_suffix);
                    cond_name = sprintf("%s-Cond%s.mat", sbj_name, save_suffix);
                    if noise(1) == 0 || noise(2) == 0
                        save(fullfile(volt_save_path, volt_name), "Umeas_NoNoise", "Uall_NoNoise", "cur_pat", "perim_mm", "noise", "flags", "volt_metadata", "-v7.3")
                    else
                        save(fullfile(volt_save_path,volt_name), "Umeas", "Uall", "cur_pat", "perim_mm", "noise", "flags", "volt_metadata", "-v7.3")
                    end
                    save(fullfile(cond_save_path,cond_name), "sigma", "sigma_GT", "nodes", "E_connect", "flags", "cond_metadata", "-v7.3")
                
        % ----------------------------------------------------------------------- %
        %%                                Plotting                                %
        % ----------------------------------------------------------------------- %
                    % Plot voltages
                    if flags.plot_volts > 0
                        figure()
                        ax1 = subplot(1,2,1);
                            p1 = scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 10, real(Uall(1:size(nodes,1),flags.plot_volts)), 'filled');
                            title(sprintf("Voltage: CP %d", flags.plot_volts))
                            xlabel("x (mm)")
                            ylabel("y (mm)")
                            zlabel("z (mm)")
                            c = colorbar;
                            c.Label.String = "mV";
                            c.Location     = 'southoutside';
                    
                        ax2 = subplot(1,2,2);
                            scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 10, imag(Uall(1:size(nodes,1),flags.plot_volts)))
                            title(sprintf("Imaginary Voltage: CP %d", flags.plot_volts))
                            xlabel("x (mm)")
                            ylabel("y (mm)")
                            zlabel("z (mm)")
                            c = colorbar;
                            c.Label.String = "mV";
                            c.Location     = 'southoutside';
                    
                        % Add a slider for x-axis min limit
                        uicontrol('Style', 'slider', 'Min', min(nodes(:,1)), 'Max', max(nodes(:,1))-1, 'Value', 0, ...
                                  'Position', [210, 30, 150, 20], ...
                                  'Callback', @(src, event) update_xlim([ax1, ax2], src.Value, max(nodes(:,1))));
                
                        % Add a slider for y-axis min limit
                        uicontrol('Style', 'slider', 'Min', min(nodes(:,2)), 'Max', max(nodes(:,2))-1, 'Value', 0, ...
                                  'Position', [210, 10, 150, 20], ...
                                  'Callback', @(src, event) update_ylim([ax1, ax2], src.Value, max(nodes(:,1))));
                    
                        if flags.do_pauses == 1
                            fprintf("      Press any key if correct\n")
                            pause
                        end
                    end % end plot_volts
                end % end running a solve
            end % end looping over heart frames
        end % end looping over permutations
    end % end looping through conditions
end % end function as a whole
    
% ----------------------------------------------------------------------- %
%%                            Custom Functions                            %
% ----------------------------------------------------------------------- %

function update_xlim(axes, minVal, maxVal)
    % Set new x-axis limits based on slider value
    for ax = axes
        ax.XLim = [minVal, maxVal];
    end
end

function update_ylim(axes, minVal, maxVal)
    % Set new y-axis limits based on slider value
    for ax = axes
        ax.YLim = [minVal, maxVal];
    end
end
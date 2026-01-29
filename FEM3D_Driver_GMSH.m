%{
Run a 3D FEM simulation on a single subject
The driver expects the subject to have the origin at the bottom, posterior,
left corner of the torso. Positive y is towards the anterior side.
All units NEED to be in standard units to work.
10/7/24 - Kyler Howard

load: tetmesh - Mesh structure containing the nodes, connectivity, and labels
load: current_pat - Current patterns for belt/patch electodes
load: sbj_info - Structure containing heights and other info for each subject

save: Umeas - Measured voltages on each electrode for each current pattern
%}

close all
clearvars -except msh_path
clc
pause('on')
start_time = tic;

% Add paths in a way that works for macs as well
addpath(fullfile("Modified_OOEIT", "ForwardProblemSolvers"))
addpath(fullfile("Modified_OOEIT", "MiscClasses"))
addpath(fullfile("Utility_Functions"))

% msh_path = "Meshes";
% msh_name = "R1053_Mesh_NoBones.mat"; % 4 months
% msh_name = "R1053_Mesh_NoBones_GMSH.mat"; % 4 months
% msh_name = "R1035_Mesh_Cylinder.mat"; % Adjust labels lines 103-112
% msh_name = "R1035_Mesh_NoBones.mat"; % 13 months
% msh_name = "R1035_Mesh_NoBones_Fine.mat"; % 13 months
% msh_name = "R1035_Mesh.mat"; % 13 months
% msh_name = "R1133_Mesh_NoBones.mat"; % 16 years, KYLER SIZED
% msh_name = "R1002_Mesh_NoBones.mat";
% msh_name = "R1043_Mesh_NoBones.mat";
% msh_name = "case101819_Mesh_Trimmed.mat";
% msh_name = "case101819_Mesh_Trimmed_NoBones.mat";
% msh_name = "case101819_Mesh_NoBones_Eroded_Rotated.mat";
% msh_name = "R1044_Mesh_NoBones_Eroded_GMSH10.mat";


if exist("msh_path", "var") && ischar(msh_path)
    [msh_name, msh_path] = uigetfile( msh_path, "Select Mesh File");
else
    [msh_name, msh_path] = uigetfile("Select Mesh File");
end

% Extracting the subject name
parts    = split(msh_name, '_');
sbj_name = parts{1};
fprintf("Running %s\n", sbj_name)

% Load the mesh
tetmesh     = load(fullfile(msh_path,msh_name), "tetmesh").tetmesh;

% ----------------------------------------------------------------------- %
%%                                Settings                                %
% ----------------------------------------------------------------------- %
% User settings
flags.do_pauses       = 0; % Decide to include pauses to check things or not
flags.solve_problem   = 1; % Decide if you want to setup (0), or fully solve (1)
flags.use_GE          = 1; % Decide if you want to use GE (1) or ACT5 (0) current patterns/conductivities
flags.do_parfor       = 1; % Decide if you want to paralize (1) or not (0)
flags.inject_current  = 0; % Decide if you want to inject ANY current (1) or only measure voltages (0)
flags.heart_BCs       = 1; % Decide if you want to include heart BCs (1) or not (0)
flags.save_heart_mesh = 0; % Decide if you want to generate and save a heart mesh (1) or not (0)
flags.do_beeps        = 0; % Decide if you want the code to beep after each simulation (1) or not (0)
flags.verbose         = 1; % Decide if you want to print status updates along the way (1) or not (0)

% Conductivity Settings
flags.set_complex     = 1; % Choice of complex (1) or real (0) conductivities
flags.const_body      = 0; % Decide if you want a solid/constant body (1) or not (0)
flags.esoph_intubate  = 0; % Decide if the esophagus is intubated (1) or not (0)
flags.max_inspiration = 0.5; % Decide if the lungs should be at inspiration (1), expiration (0), or somewhere in-between
flags.equal_vent      = 1; % Decide if you want equal ventilation (1) or split (0)
flags.left_only       = 0; % Decide if you want only ventilation on the left side (1) or not (0)
flags.right_only      = 0; % Decide if you want only ventilation on the right side (1) or not (0)
flags.permute_conds   = 0; % Decide if you want random conds (1) or not (0)

% Plot settings
flags.plot_slices     = 0; % Plot individual slices when going slice by slice
flags.plot_trachea    = 1; % Plotting of carina height & trachea orientation
flags.plot_electrodes = 1; % Plotting of electrode consturction
flags.plot_conds      = 0; % Plotting of conductivities
flags.plot_GTs        = 1; % Plot ground truth images
flags.plot_internal   = 0; % Plotting of internal nodes
flags.plot_volts      = 1; % Plotting of nodal voltages
flags.fixed_range     = 1; % Set GT plots to be a standard range

flags.CP_choice       = 1; % Choice of current pattern for patches
    % 1: Standard pattern
    % 2: 4x8 pattern
flags.E_choice        = 4; % Choice of Electrode configuration
    % 1: Large patch front back  (GE Patch)
    % 2: Small patch front back  (GE Patch)
    % 3: Two rows of large belts (GE Belt)
    % 4: Two rows of small belts (GE Belt)
    % 5: Custom electrodes

% Custom Electrode Settings
flags.E_type          = "belt";   % Choice between "patch" and "belt"
flags.E_shape         = "circle"; % Choice between "circle" and "rectangle"
flags.E_dia           = 30;  % Diameter of electrode in mm (for circle)
flags.E_width         = 22;  % Width  of electrode in mm (for rectangle)
flags.E_height        = 29;  % Height of electrode in mm (for rectangle)
flags.gap_width       = 46.675; % Gap between electrodes horizontally in mm (edge-edge) (for patch) %2.5 / 46.675
flags.gap_height      = 32.875; % Gap between electrodes vertically in mm (edge-edge) (for patch) %2.5 / 32.3875
flags.E_count         = [16];  % Number of electrodes per row (for belt), or matrix of how many rows and columns (for patch)

% Save the correct electrode settings, not the custom ones when using the
% standard settings
flags = Construct_Electrode_Settings(flags);

noise = [0, 0, 0, 0];       % Noise and error parameters
    % noise_rel = err(1);
    % noise_abs = err(2);
    % e_systematic_rel = err(3);
    % e_systematic_abs = err(4);

if contains(msh_name, 'NoBones')
    flags.are_bones = 0; % There are no bones
else
    flags.are_bones = 1; % There are bones
end

if isempty(gcp('nocreate')) && flags.solve_problem == 1 && flags.do_parfor == 1
    % Open a parallel pool
    parpool;
end
% ----------------------------------------------------------------------- %
%%                                 Setup                                  %
% ----------------------------------------------------------------------- %
if flags.verbose == 1
    fprintf("   Defining Variables\n")
end
% Splitting the mesh structure
connectivity = double(tetmesh.cell'); % indices of nodes making up an element
nodes        = double(tetmesh.node'); % xyz coordinates in mm
labels       = tetmesh.field';        % material id for each element
% clear tetmesh parts

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
    % flags.lungmap_set = 1;
elseif contains(sbj_name, "case1")
    sbj_sheet = readtable("CT Data Boundaries.xlsx", "Sheet","New Mexico Baby Set");
    % flags.lungmap_set = 0;
elseif contains(sbj_name, "EIT1")
    sbj_sheet = readtable("CT Data Boundaries.xlsx", "Sheet","Anschutz Set");
    % flags.lungmap_set = 0;
elseif contains(sbj_name, "MCR0")
    sbj_sheet = readtable("CT Data Boundaries.xlsx", "Sheet","MCR Set");
    % flags.lungmap_set = 0;
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

% DELETE ME!!! - Manually erasing spinal cords for those who have it while
% testing
if size(unique(labels),1) > external
    labels(labels == external - 1) = soft_tissue - 1;
    labels(labels == external)     = external - 1;
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
[nodes, sbj_info] = Rotate_and_Translate_Body_GMSH(nodes, organ_connects, carina_height, T5_height, T8_height, flags);

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
% carina_height = ((max(boundary_nodes(:,3)) + min(boundary_nodes(:,3)))/2) + 15/2 + 5/2; % DELETE ME!!!! FOR CYLINDER ONLY
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

num_E_faces = zeros(1,L);
E_areas     = zeros(1,L);
for i = 1:length(E_connect)
    num_E_faces(i) = size(E_connect{i},1);
    
    for j = 1:num_E_faces(i)
        face_nodes = E_connect{i}(j,:);
        point1 = nodes(face_nodes(1),:);
        point2 = nodes(face_nodes(2),:);
        point3 = nodes(face_nodes(3),:);

        E_areas(i) = E_areas(i) + norm(cross(point3-point1, point3-point2)) / 2;
    end
end
clear i j

if flags.E_choice < 3 || (flags.E_choice == 5 && flags.E_type == "patch")
    front = 1:16;
    back  = 17:32;
else
    front = [5:12, 21:28];
    back  = [1:4,13:20,29:32];
end

if flags.plot_electrodes == 1
    figure()
        subplot(1,2,1)
            hold on
            bar(front,num_E_faces(front),'EdgeColor',"#0072BD",'FaceColor',"#0072BD")
            bar(back,num_E_faces(back),  'EdgeColor',"#4DBEEE",'FaceColor',"#4DBEEE")
            plot(1:L, mean(num_E_faces)*ones(1,L),'r')
            xlabel("Electrode")
            title("Number of Faces")
            title(sprintf("Number of Faces\nMean: %.0f, STD: %.0f", mean(num_E_faces), std(num_E_faces)))
            legend("Front", "Back")
        subplot(1,2,2)
            hold on
            bar(front,E_areas(front),'EdgeColor',"#0072BD",'FaceColor',"#0072BD")
            bar(back,E_areas(back),  'EdgeColor',"#4DBEEE",'FaceColor',"#4DBEEE")
            plot(1:L, mean(E_areas)*ones(1,L),'r')
            xlabel("Electrode")
            title(sprintf("Electrode Area (mm²)\nMean: %.2f, STD: %.2f", mean(E_areas), std(E_areas)))
            legend("Front", "Back")
end

carina_plane = find(boundary_nodes(:,3) < sbj_info.carina + 1.25 & boundary_nodes(:,3) > sbj_info.carina - 1.25);

if flags.plot_electrodes == 1
    figure()
        hold on
        scatter3(boundary_nodes(:,1), boundary_nodes(:,2), boundary_nodes(:,3), 'MarkerEdgeColor', [0.3010 0.7450 0.9330], 'LineWidth',0.05)
        % for cell_i = 1:length(E_nodes)
        %     scatter3(E_nodes{cell_i}(:,1), E_nodes{cell_i}(:,2), E_nodes{cell_i}(:,3), 'r')
        % end
        % scatter3(boundary_nodes(carina_plane,1), boundary_nodes(carina_plane, 2), boundary_nodes(carina_plane,3),'b','filled')
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
% return

% ----------------------------------------------------------------------- %
%%                            Heart Conditions                            %
% ----------------------------------------------------------------------- %
% Setting up heart boundary condition information
if flags.heart_BCs == 1
    heart_BCs = load(fullfile("Heart_BCs", sprintf("%s_Heart_BC_frame_1.mat",sbj_name)));
    
    heart_BC_active_current = heart_BCs.current_density(1)*10^-3; % CONVERT mAMPS TO AMPS
    heart_BC_middle_current = heart_BCs.current_density(2)*10^-3; % CONVERT mAMPS TO AMPS
    heart_BC_active_indices = heart_BCs.nodeGroups.active;
    heart_BC_middle_indices = heart_BCs.nodeGroups.middle;
    clear heart_BCs

    % Convert the given indices into the node coordinates
    heart_surface_nodes   = load(fullfile("Heart_Meshes", sprintf("%s_Heart_Mesh.mat",sbj_name))).heart_surface_nodes;
    heart_faces           = load(fullfile("Heart_Meshes", sprintf("%s_Heart_Mesh.mat",sbj_name))).heart_faces;
    heart_BC_active_nodes = {heart_surface_nodes(heart_BC_active_indices,:)};
    heart_BC_middle_nodes = {heart_surface_nodes(heart_BC_middle_indices,:)};

    % Find the heart mesh faces that make up the BCs
    heart_active_connect = Align_Electrode_Faces(nodes, heart_faces, heart_BC_active_nodes, flags);
    heart_middle_connect = Align_Electrode_Faces(nodes, heart_faces, heart_BC_middle_nodes, flags);

    % Append everything to the electrode BCs
    E_connect = [E_connect, heart_active_connect, heart_middle_connect];
    cur_pat   = [cur_pat; heart_BC_active_current*ones(length(heart_active_connect), K); heart_BC_middle_current*ones(length(heart_middle_connect), K)]; 

    % Update L/K to reflect the new BC
    L = L + length(heart_active_connect) + length(heart_middle_connect);
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
    fprintf("      It took %.2f minutes to initialize solver\n", init_time/60)
            
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
    if flags.heart_BCs == 1
        flags.zeta(end-length(heart_active_connect)-length(heart_middle_connect)+1:end) = 0.116;
    end
    solver.mode = 'current';
    solver.zeta = flags.zeta;
end

% ----------------------------------------------------------------------- %
%%                         Assign Conductivities                          %
% ----------------------------------------------------------------------- %
% for i = 1%:6
    % if i == 1
    %     flags.esoph_intubate  = 0; % Decide if the esophagus is intubated (1) or not (0)
    %     flags.max_inspiration = 0.5; % Decide if the lungs should be at inspiration (1), expiration (0), or in-between (0.5)
    %     flags.left_only       = 0; % Decide if you want only ventilation on the left side (1) or not (0)
    %     flags.right_only      = 0; % Decide if you want only ventilation on the right side (1) or not (0)
    % elseif i == 2
    %     flags.esoph_intubate  = 0; % Decide if the esophagus is intubated (1) or not (0)
    %     flags.max_inspiration = 0.2; % Decide if the lungs should be at inspiration (1), expiration (0), or in-between (0.5)
    %     flags.left_only       = 0; % Decide if you want only ventilation on the left side (1) or not (0)
    %     flags.right_only      = 0; % Decide if you want only ventilation on the right side (1) or not (0)
    % elseif i == 3
    %     flags.esoph_intubate  = 0; % Decide if the esophagus is intubated (1) or not (0)
    %     flags.max_inspiration = 0; % Decide if the lungs should be at inspiration (1), expiration (0), or in-between (0.5)
    %     flags.left_only       = 0; % Decide if you want only ventilation on the left side (1) or not (0)
    %     flags.right_only      = 0; % Decide if you want only ventilation on the right side (1) or not (0)
    % elseif i == 4
    %     flags.esoph_intubate  = 0; % Decide if the esophagus is intubated (1) or not (0)
    %     flags.max_inspiration = 1; % Decide if the lungs should be at inspiration (1), expiration (0), or in-between (0.5)
    %     flags.left_only       = 1; % Decide if you want only ventilation on the left side (1) or not (0)
    %     flags.right_only      = 0; % Decide if you want only ventilation on the right side (1) or not (0)
    % elseif i == 5
    %     flags.esoph_intubate  = 0; % Decide if the esophagus is intubated (1) or not (0)
    %     flags.max_inspiration = 1; % Decide if the lungs should be at inspiration (1), expiration (0), or in-between (0.5)
    %     flags.left_only       = 0; % Decide if you want only ventilation on the left side (1) or not (0)
    %     flags.right_only      = 1; % Decide if you want only ventilation on the right side (1) or not (0)
    % end

    if flags.verbose == 1
        fprintf("   Assigning Conductivities\n")
    end

    % Creating the conductivity vector at the nodes (IN SIEMENS PER METER)
    sigma = Assign_Conductivities(nodes, connectivity, labels, lung_nodes, flags);

% ----------------------------------------------------------------------- %
%%                       Create Ground Truth Images                       %
% ----------------------------------------------------------------------- %
    t_GT = 5;
    E_heights = cellfun(@(x) mean(x(:,3)), E_nodes);
    if flags.E_type == "patch"
        if flags.CP_choice == 1
            E_heights = sort([mean(E_heights([1:4,17:20])), mean(E_heights([5:8,21:24])), mean(E_heights([9:12,25:28])), mean(E_heights([13:16,29:32]))], "descend");
        else
            E_heights = sort([mean(E_heights(1:8)), mean(E_heights(9:16)), mean(E_heights(17:24)), mean(E_heights(25:32))], "descend");
        end
    else
        E_heights = sort([mean(E_heights(1:16)), mean(E_heights(17:32))], "descend");
    end
    
    if flags.E_type == "patch"
        sigma_GT = cell(4,1);
    else
        sigma_GT = cell(2,1);
    end
    for row = 1:length(E_heights)
        plane = (nodes(:,3) > E_heights(row) - t_GT/2) & (nodes(:,3) < E_heights(row) + t_GT/2);
        plane = plane & sigma~=0;
    
        sigma_GT{row} = plane;
    end
        
    if flags.plot_GTs == 1
        if flags.set_complex == 1
            num_col = 2;
        else
            num_col = 1;
        end

        % Initialize
        shift = 0;
        rowStart = 1;
        plane_nodes = zeros(sum(cellfun(@sum, sigma_GT)),2);
        plane_sigma = zeros(sum(cellfun(@sum, sigma_GT)),1);
        for k = 1:length(E_heights)
            % Current matrices
            nodemat = [nodes(sigma_GT{k},1), nodes(sigma_GT{k},2)];
            sigmat  = sigma(sigma_GT{k});
            nrows   = size(nodemat,1);
            PN = nodemat;
    
            % Apply shift to second column
            PN(:,2) = PN(:,2) - shift;
    
            % Append to full matrices
            plane_nodes(rowStart:rowStart+nrows-1, :) = PN;
            plane_sigma(rowStart:rowStart+nrows-1, :) = sigmat;
    
            % Update shift: add the max of y values plus 33%
            shift = shift + 4/3*max(nodemat(:,2));
            rowStart = rowStart + nrows;
        end

        figure('color','w','Position',[573,337.67,700,420]);
        for col = 1:num_col
            subplot(1,num_col,col)
            
            if col == 1
                scatter(plane_nodes(:,1), plane_nodes(:,2), 10, real(plane_sigma), 'filled')
            else
                scatter(plane_nodes(:,1), plane_nodes(:,2), 10, imag(plane_sigma), 'filled')
            end
    
            % Flip the x axis horizontally to be in DICOM standard
            % set(gca, 'XDir', 'reverse'); % KH: it already was I think
            axis equal off
    
            colormap("jet")
            if flags.fixed_range == 1
                clim([0, 0.8])
            else
                if col == 1
                    clim([min(real(plane_sigma)), max(real(plane_sigma))])
                else
                    clim([min(imag(plane_sigma)), max(imag(plane_sigma))])
                end
            end
        
            if col == 1
                title(sprintf("Ground Truth\nConductivity"))
            else
                title(sprintf("Ground Truth\nSusceptivity"))
            end

            c = colorbar("eastoutside");
            c.Label.String = "S/m";
        end

        % %%
        % figure('color','w','Position',[573,337.67,700,420]);
        % for col = 1:num_col
        %     for row = 1:length(sigma_GT)
        %         if flags.E_type == "patch"
        %             subplot(4,num_col,2*row + col - 2)
        %         else
        %             subplot(2,num_col,2*row + col - 2)
        %         end
        % 
        %         if col == 1
        %             scatter(nodes(sigma_GT{row},1), nodes(sigma_GT{row},2), 10, real(sigma(sigma_GT{row})), 'filled')
        %         else
        %             scatter(nodes(sigma_GT{row},1), nodes(sigma_GT{row},2), 10, imag(sigma(sigma_GT{row})), 'filled')
        %         end
        % 
        %         % Flip the x axis horizontally to be in DICOM standard
        %         % set(gca, 'XDir', 'reverse'); % KH: it already was I think
        %         axis equal off
        % 
        %         colormap("jet")
        %         if flags.fixed_range == 1
        %             clim([0, 0.8])
        %         else
        %             if col == 1
        %                 clim([min(real(sigma)), max(real(sigma))])
        %             else
        %                 clim([min(imag(sigma)), max(imag(sigma))])
        %             end
        %         end
        % 
        %         if row == 1
        %             if col == 1
        %                 title(sprintf("Ground Truth\nConductivity"))
        %             else
        %                 title(sprintf("Ground Truth\nSusceptivity"))
        %             end
        %         end
        % 
        %     end
        %     if flags.set_complex == 1
        %         if col == 1
        %             c = colorbar('Position', [0.48 0.11 0.025 0.815]);
        %             c.Label.String = "S/m";
        %         else
        %             c = colorbar('Position', [0.91 0.11 0.025 0.815]);
        %             c.Label.String = "S/m";
        %         end
        %     else
        %         c = colorbar('Position', [0.7 0.11 0.025 0.815]);
        %         c.Label.String = "S/m";
        %     end
        % end
    end

% ----------------------------------------------------------------------- %
%%                                 Solve                                  %
% ----------------------------------------------------------------------- %

% zeta = [0.041];
% for ii = 1:length(zeta)
%     solver.zeta = zeta(ii)*ones(L,1);
%     flags.zeta  = solver.zeta;    

    if flags.solve_problem == 1
        fprintf("   Solving Forward Problem\n")
    
        % Set the current pattern (IN AMPS)
        solver.Iel = cur_pat;
    
        % Solve the voltage
        solve_start          = tic;
        if flags.do_parfor == 1
            [Umeas, Imeas, Uall] = MF_Simulation(fmesh, [], sigma, solver.zeta, solver, 'current', noise);
        else
            [Umeas, Imeas, Uall] = MF_Simulation2(fmesh, [], sigma, solver.zeta, solver, 'current', noise);
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
        save_suffix = Make_Save_Name(solver.zeta, flags);

        % Create metadata to include with the saving
        volt_metadata = Make_Metadata("Volt");
        cond_metadata = Make_Metadata("Cond");
    
        % Save the voltages and conductivties
        volt_name = sprintf("%s_Volt%s.mat", sbj_name, save_suffix);
        cond_name = sprintf("%s_Cond%s.mat", sbj_name, save_suffix);
        if noise(1) == 0 || noise(2) == 0
            save(fullfile("Results", volt_name), "Umeas_NoNoise", "Uall_NoNoise", "cur_pat", "perim_mm", "noise", "flags", "volt_metadata", "-v7.3")
        else
            save(fullfile("Results",volt_name), "Umeas", "Uall", "cur_pat", "perim_mm", "noise", "flags", "volt_metadata", "-v7.3")
        end
        save(fullfile("Results",cond_name), "sigma", "sigma_GT", "nodes", "E_connect", "flags", "cond_metadata", "-v7.3")
    
    % ----------------------------------------------------------------------- %
    %%                                Plotting                                %
    % ----------------------------------------------------------------------- %
    
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
        end
    end
% end

stop_time = toc(start_time);
if flags.verbose == 1
    fprintf("\n   It took %.2f hours to solve the forward problem\n", stop_time / 3600)
end

if flags.do_beeps == 1
    beep
end

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
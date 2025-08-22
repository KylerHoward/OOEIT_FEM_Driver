%{
Run a 3D FEM simulation on a single subject
10/7/24 - Kyler Howard

load: tetmesh - Mesh structure containing the nodes, connectivity, and labels
load: current_pat - Current patterns for belt/patch electodes
load: sbj_info - Structure containing heights and other info for each subject

save: Umeas - Measured voltages on each electrode for each current pattern
%}

close all
clear
clc
pause('on')
start_time = tic;

addpath Modified_OOEIT\ForwardProblemSolvers\ Modified_OOEIT\MiscClasses\ Utility_Functions\

msh_path = "Meshes";
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
msh_name = "case101819_Mesh_NoBones_Eroded_Rotated.mat";

% [msh_name, msh_path] = uigetfile("Select Mesh File"); 

% Load the mesh
tetmesh     = load(fullfile(msh_path,msh_name), "tetmesh").tetmesh;

% ----------------------------------------------------------------------- %
%%                                Settings                                %
% ----------------------------------------------------------------------- %
% User settings
flags.do_pauses       = 0; % Decide to include pauses to check things or not
flags.solve_problem   = 0; % Decide if you want to setup (0), or fully solve (1)
flags.use_GE          = 1; % Decide if you want to use GE (1) or ACT5 (0) current patterns/conductivities
flags.do_parfor       = 0; % Decide if you want to paralize (1) or not (0)

% Conductivity Settings
flags.set_complex     = 0; % Choice of complex (1) or real (0) conductivities
flags.const_body      = 0; % Decide if you want a solid/constant body (1) or not (0)
flags.esoph_intubate  = 0; % Decide if the esophagus is intubated (1) or not (0)
flags.max_inspiration = 1; % Decide if the lungs should be at inspiration (1), expiration (0), or somewhere in-between
flags.equal_vent      = 1; % Decide if you want equal ventilation (1) or split (0)
flags.left_only       = 0; % Decide if you want only ventilation on the left side (1) or not (0)
flags.right_only      = 0; % Decide if you want only ventilation on the right side (1) or not (0)
flags.permute_conds   = 0; % Decide if you want random conds (1) or not (0)

% Plot settings
flags.plot_slices     = 0; % Plot individual slices when going slice by slice
flags.plot_trachea    = 1; % Plotting of carina height
flags.plot_electrodes = 1; % Plotting of electrode consturction
flags.plot_conds      = 1; % Plotting of conductivities
flags.plot_internal   = 0; % Plotting of internal nodes
flags.plot_volts      = 1; % Plotting of nodal voltages

flags.CP_choice       = 1; % Choice of current pattern
    % 1: Standard pattern
    % 2: 4x8 pattern
flags.E_choice        = 3; % Choice of Electrode configuration
    % 1: Large patch front back  (GE Patch)
    % 2: Small patch front back  (GE Patch)
    % 3: Two rows of large belts (GE Belt)
    % 4: Two rows of small belts (GE Belt)
    % 5: Custom electrodes

% Custom Electrode Settings
flags.E_type          = "patch";   % Choice between "patch" and "belt"
flags.E_shape         = "circle"; % Choice between "circle" and "rectangle"
flags.E_dia           = 20;  % Diameter of electrode in mm (for circle)
flags.E_width         = 22;  % Width  of electrode in mm (for rectangle)
flags.E_height        = 29;  % Height of electrode in mm (for rectangle)
flags.gap_width       = 46.675; % Gap between electrodes vertically in mm (edge-edge) (for patch) %2.5 / 46.675
flags.gap_height      = 32.875; % Gap between electrodes vertically in mm (edge-edge) (for patch) %2.5 / 32.3875
flags.E_count         = [4,4];  % Number of electrodes per row (for belt), or matrix of how many rows and columns (for patch)

err = [0, 0, 0, 0];       % Noise and error parameters
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
    parpool
end
% ----------------------------------------------------------------------- %
%%                                 Setup                                  %
% ----------------------------------------------------------------------- %
fprintf("Defining Variables\n")
% Splitting the mesh structure
connectivity = double(tetmesh.cell'); % indices of nodes making up an element
nodes        = double(tetmesh.node'); % xyz coordinates in mm
labels       = tetmesh.field';        % material id for each element

if (flags.E_choice == 5 && flags.E_type == "patch") || flags.E_choice <= 2
    if flags.CP_choice == 1
        current_pattern = load(fullfile("Current_Patterns", "Box4by4aMxStrength8.mat"), "Box4by4aMxStrength8").Box4by4aMxStrength8;
    elseif flags.CP_choice == 2
        current_pattern = load(fullfile("Current_Patterns", "CP32_4x8.mat"), "CP32_4x8").CP32_4x8;
    end
    L               = size(current_pattern, 1); % Number of electrodes
    num_CP          = size(current_pattern, 2); % Number of current patterns
elseif (flags.E_choice == 5 && flags.E_type == "belt") || flags.E_choice > 2
    % current_pat = load(fullfile("Current_Patterns", "CP32_16x2_M1.mat"), "Cur_pat3D").Cur_pat3D; % Normalized CP
    if flags.use_GE == 1
        % Current pattern from GE is in µA
        current_pattern = load(fullfile("Current_Patterns", "clinical_belt_CP.mat"), "CP_belt_16x2").CP_belt_16x2;
    else
        % Current pattern from ACT5 is in ???
        error("ACT5 Current Not Loaded")
    end

    L      = size(current_pattern, 1); % Number of electrodes
    num_CP = size(current_pattern, 2); % Number of current patterns

    % Scale the current pattern
    CP_scale        = load(fullfile("Current_Patterns", "clinical_belt_CP.mat"), "CPscale_belt_16x2").CPscale_belt_16x2;
    CP_scale        = repmat(CP_scale', 1, num_CP);
    current_pattern = CP_scale .* current_pattern;
end

% Extracting the subject name
parts    = split(msh_name, '_');
sbj_name = parts(1);
clear tetmesh parts

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
upside_down   = sbj_sheet.UpsideDown(row_num);
facing_dir    = sbj_sheet.PositiveYIs(row_num);
carina_height = sbj_sheet.CarinaFromBottom_mm(row_num);
T5_height     = sbj_sheet.T5FromBottom_mm(row_num);
T8_height     = sbj_sheet.T8FromBottom_mm(row_num);

% %KH DELETE ME TESTING FOR R1133
if contains(msh_name, "R1133")
    carina_height = carina_height - 50;
end
if contains(msh_name, "case101819")
    carina_height = 150;
end

sbj_info.carina = carina_height;
sbj_info.T5     = T5_height;
sbj_info.T8     = T8_height;

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

% ----------------------------------------------------------------------- %
%%                        Rotation and Translation                        %
% ----------------------------------------------------------------------- %
% Determine if the body is upside down or on it's side
[body_nodes, ~]     = Get_Tet_Nodes(nodes, organ_connects{soft_tissue});
[trachea_nodes, ~]  = Get_Tet_Nodes(nodes, organ_connects{trachea});

% Check if the body is tilted on it's side
if (median(trachea_nodes(:,2)) - mean(trachea_nodes(:,2))) > 0.05*range(trachea_nodes(:,2))
    fprintf("Sitting them upright\n")
    theta = -pi/2;
    rotationMatrix = [1, 0,          0;...
                      0, cos(theta), sin(theta);...
                      0, sin(theta), cos(theta)];

    % Rotate and shift the body
    nodes      = nodes*rotationMatrix;
    nodes(:,1) = nodes(:,1) * -1;
    nodes(:,3) = nodes(:,3) + abs(min(nodes(:,3)));

    % Recalculate shifted nodes
    [body_nodes, ~]    = Get_Tet_Nodes(nodes, organ_connects{soft_tissue});
    [trachea_nodes, ~] = Get_Tet_Nodes(nodes, organ_connects{trachea});
end

% Check the excel doc
% if upside_down == 0
% Check if the difference between the trachea and the top of the body is less than 5% of the total height
if abs(min(trachea_nodes(:,3)) - min(body_nodes(:,3))) < 0.05*max(nodes(:,3))
    fprintf("Rotating Body Upright\n")
    theta = pi;
    rotationMatrix = [cos(theta), 0, -sin(theta);...
                      0,          1,  0;...
                      sin(theta), 0,  cos(theta)];

    % Rotate and shift the body
    nodes      = nodes*rotationMatrix;
    nodes(:,1) = nodes(:,1) * -1;
    nodes(:,3) = nodes(:,3) + abs(min(nodes(:,3)));

    % Recalculate shifted nodes
    [body_nodes, ~]    = Get_Tet_Nodes(nodes, organ_connects{soft_tissue});
    [trachea_nodes, ~] = Get_Tet_Nodes(nodes, organ_connects{trachea});
end

[heart_nodes, ~]   = Get_Tet_Nodes(nodes, organ_connects{heart});
heart_center       = range(heart_nodes(:,2))/2 + min(heart_nodes(:,2));
body_center        = range(body_nodes(:,2))/2 + min(body_nodes(:,2));
trachea_center     = range(trachea_nodes(:,2))/2 + min(trachea_nodes(:,2));

% if facing_dir == "Posterior"
% if (heart_center < body_center) || (trachea_center > body_center) % Check if the heart/trachea is on the wrong side
if range(body_nodes(:,1)) < range(body_nodes(:,2)) % Check if the heart/trachea is on the wrong side
    fprintf("Rotating the Chest Forward\n")
    % theta = pi;
    theta = -pi/2;
    rotationMatrix = [cos(theta), -sin(theta), 0;...
                      sin(theta),  cos(theta), 0;...
                      0,           0,          1];

    % Rotate and shift the body
    nodes      = nodes*rotationMatrix;
    nodes(:,1) = nodes(:,1) + abs(min(nodes(:,1)));
    nodes(:,2) = nodes(:,2) + abs(min(nodes(:,2)));
end

fprintf("Translating to the Origin\n")
[body_nodes, ~] = Get_Tet_Nodes(nodes, organ_connects{soft_tissue});
nodes           = nodes         - min(body_nodes);
sbj_info.carina = carina_height - min(body_nodes(:,3));
sbj_info.T5     = T5_height     - min(body_nodes(:,3));
sbj_info.T8     = T8_height     - min(body_nodes(:,3));
[body_nodes, ~] = Get_Tet_Nodes(nodes, organ_connects{soft_tissue});

if contains(sbj_name, "case101819")
    sbj_info.carina = 150;
    sbj_info.T5 = 165;
    sbj_info.T8 = 100;
end

% ----------------------------------------------------------------------- %
%%                              Surface Nodes                             %
% ----------------------------------------------------------------------- %
fprintf("Extracting Trachea and Surface Nodes\n")
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

clear trachea_faces heart_faces esophagus_faces
clear lung_intersect trachea_intersect heart_intersect esophagus_intersect inner_intersects surface_faces

if flags.plot_trachea == 1
    figure()
    subplot(1,2,1)
        scatter3(trachea_nodes(:,1), trachea_nodes(:,2), trachea_nodes(:,3), "MarkerEdgeAlpha", 0.2)
        hold on
        scatter3(mean(trachea_nodes(:,1)), mean(trachea_nodes(:,2)), sbj_info.carina, 'r', 'filled')
        scatter3(trachea_nodes(startsWith(string(trachea_nodes(:,3)), sprintf("%.1f", sbj_info.carina)), 1), trachea_nodes(startsWith(string(trachea_nodes(:,3)), sprintf("%.1f", sbj_info.carina)), 2), sbj_info.carina, 'r', 'filled')
        % constantplane('z', sbj_info.carina) R2024b
        title(sprintf("Carina: %.2f mm", sbj_info.carina))
    subplot(1,2,2)
        pdeplot3D(nodes', organ_connects{lung}')
end

% ----------------------------------------------------------------------- %
%%                             Make Electrodes                            %
% ----------------------------------------------------------------------- %
fprintf("Making Electrodes\n")
tic
% carina_height = ((max(boundary_nodes(:,3)) + min(boundary_nodes(:,3)))/2) + 15/2 + 5/2; % DELETE ME!!!! FOR CYLINDER ONLY
[E_nodes, perim_mm] = Make_Electrodes3(boundary_nodes, nodes, body_faces, sbj_info, flags);
if iscell(E_nodes) == 0
    E_nodes = {E_nodes};
end

% Find the surface mesh faces that make up the electrodes
E_connect = Align_Electrode_Faces(nodes, body_faces, E_nodes, flags);

make_time = toc;
fprintf("   It took %.2f seconds to make the electrodes\n", make_time)

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
            legend("Front", "Back")
        subplot(1,2,2)
            hold on
            bar(front,E_areas(front),'EdgeColor',"#0072BD",'FaceColor',"#0072BD")
            bar(back,E_areas(back),  'EdgeColor',"#4DBEEE",'FaceColor',"#4DBEEE")
            plot(1:L, mean(E_areas)*ones(1,L),'r')
            xlabel("Electrode")
            title(sprintf("Electrode Area (mm²), STD: %.2f mm²", std(E_areas)))
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
        for i = 1:length(E_connect)
            trimesh(E_connect{i}, nodes(:,1), nodes(:,2), nodes(:,3), 'FaceColor', 'r', 'edgeColor', 'r');
        end
        xlabel('X (mm)')
        ylabel('Y (mm)')
        zlabel('Z (mm)')
        axis equal
    if flags.do_pauses == 1
        fprintf("   Press any key if correct\n")
        pause
    end
end
return
% ----------------------------------------------------------------------- %
%%                             Create Solver                              %
% ----------------------------------------------------------------------- %
% for i = 1:5
%     if i == 1
%         flags.esoph_intubate  = 0; % Decide if the esophagus is intubated (1) or not (0)
%         flags.max_inspiration = 0.5; % Decide if the lungs should be at inspiration (1), expiration (0), or in-between (0.5)
%         flags.left_only       = 0; % Decide if you want only ventilation on the left side (1) or not (0)
%         flags.right_only      = 0; % Decide if you want only ventilation on the right side (1) or not (0)
%     elseif i == 2
%         flags.esoph_intubate  = 0; % Decide if the esophagus is intubated (1) or not (0)
%         flags.max_inspiration = 1; % Decide if the lungs should be at inspiration (1), expiration (0), or in-between (0.5)
%         flags.left_only       = 0; % Decide if you want only ventilation on the left side (1) or not (0)
%         flags.right_only      = 0; % Decide if you want only ventilation on the right side (1) or not (0)
%     elseif i == 3
%         flags.esoph_intubate  = 0; % Decide if the esophagus is intubated (1) or not (0)
%         flags.max_inspiration = 0; % Decide if the lungs should be at inspiration (1), expiration (0), or in-between (0.5)
%         flags.left_only       = 0; % Decide if you want only ventilation on the left side (1) or not (0)
%         flags.right_only      = 0; % Decide if you want only ventilation on the right side (1) or not (0)
%     elseif i == 4
%         flags.esoph_intubate  = 0; % Decide if the esophagus is intubated (1) or not (0)
%         flags.max_inspiration = 1; % Decide if the lungs should be at inspiration (1), expiration (0), or in-between (0.5)
%         flags.left_only       = 1; % Decide if you want only ventilation on the left side (1) or not (0)
%         flags.right_only      = 0; % Decide if you want only ventilation on the right side (1) or not (0)
%     elseif i == 5
%         flags.esoph_intubate  = 0; % Decide if the esophagus is intubated (1) or not (0)
%         flags.max_inspiration = 1; % Decide if the lungs should be at inspiration (1), expiration (0), or in-between (0.5)
%         flags.left_only       = 0; % Decide if you want only ventilation on the left side (1) or not (0)
%         flags.right_only      = 1; % Decide if you want only ventilation on the right side (1) or not (0)
%     end

    if flags.solve_problem == 1
        % Create a mesh-object from nodes, connections, and electrode connections (NODES IN METERS)
        fprintf("Making Forward Mesh\n")
        if flags.do_parfor == 1
            fmesh = ParForwardMesh1st(nodes*1e-3, connectivity, E_connect);
        else
            fmesh = ForwardMesh1st(nodes*1e-3, connectivity, E_connect);
        end
        
        % Initialize the forward problem solver (FPS). This is the object, that
        % (utilizing the mesh-object given to it) computes forward problem solutions
        % using the FEM approximation of the complete electrode model (CEM).
        fprintf("Initializing FEM Solver\n")
        init_start = tic;
        if flags.do_parfor == 1
            solver = MF_EITFEM(fmesh);
        else
            solver = EITFEM(fmesh);
        end
        init_time  = toc(init_start);
        fprintf("   It took %.2f minutes to initialize solver\n", init_time/60)
                
        % Set solving mode to injecting current and measuring voltages
        % Set the contact impedance - had 2.4 originally in 2D FEM, then 4.8, then 0.05
        % Fine testing with closly matched rect. electrodes:  0.2475 for Subj011
        % Fine testing with closly matched patch. electrodes: 0.05 for Kyler
        solver.mode = 'current';
        solver.zeta = 0.05*ones(L,1);
    end

% ----------------------------------------------------------------------- %
%%                         Assign Conductivities                          %
% ----------------------------------------------------------------------- %
    fprintf("Assigning Conductivities\n")

    % Creating the conductivity vector at the nodes (IN SIEMENS PER METER)
    sigma = Assign_Conductivities(nodes, connectivity, labels, lung_nodes, flags);

% ----------------------------------------------------------------------- %
%%                                 Solve                                  %
% ----------------------------------------------------------------------- %
% zeta = [0.05, 0.10, 0.15];
% for i = 1:3
%     solver.zeta = zeta(i)*ones(L,1);
    
    if flags.solve_problem == 1
        fprintf("Solving Forward Problem\n")
    
        % Set the current pattern (IN AMPS)
        solver.Iel = current_pattern * 1e-6;
    
        % Solve the voltage
        solve_start          = tic;
        if flags.do_parfor == 1
            [Umeas, Imeas, Uall] = MF_Simulation(fmesh, [], sigma, solver.zeta, solver, 'current', err);
        else
            [Umeas, Imeas, Uall] = MF_Simulation2(fmesh, [], sigma, solver.zeta, solver, 'current', err);
        end
        solve_time           = toc(solve_start);
        fprintf("   It took %.2f minutes to solve for voltages\n", solve_time/60)
    
        Umeas = Umeas * 1e3; % Convert V to mV
        Uall  = Uall  * 1e3; % Convert V to mV
        Umeas = reshape(Umeas, L, num_CP);
    
        % Create the suffix for saving the file based on the settings chosen
        save_suffix = Make_Save_Name(solver.zeta, flags);
    
        % Save the voltages and conductivties
        volt_name = sprintf("Volt_%s%s.mat", sbj_name, save_suffix);
        cond_name = sprintf("Cond_%s%s.mat", sbj_name, save_suffix);
        save(fullfile("Results",volt_name), "Umeas", "Uall", "current_pattern", "perim_mm", "-v7.3")
        save(fullfile("Results",cond_name), "sigma", "nodes", "-v7.3")
    
    % ----------------------------------------------------------------------- %
    %%                                Plotting                                %
    % ----------------------------------------------------------------------- %
    
        if flags.plot_volts > 0
            figure()
            ax1 = subplot(1,2,1);
                p1 = scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 10, real(Uall(1:end-num_CP,flags.plot_volts)), 'filled');
                title(sprintf("Voltage: CP %d", flags.plot_volts))
                xlabel("x (mm)")
                ylabel("y (mm)")
                zlabel("z (mm)")
                c = colorbar;
                c.Label.String = "mV";
                c.Location     = 'southoutside';
        
            ax2 = subplot(1,2,2);
                scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 10, imag(Uall(1:end-num_CP,flags.plot_volts)))
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
                fprintf("   Press any key if correct\n")
                pause
            end
        end
    end
% end
stop_time = toc(start_time);
fprintf("\nIt took %.2f hours to solve the forward problem\n", stop_time / 3600)


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
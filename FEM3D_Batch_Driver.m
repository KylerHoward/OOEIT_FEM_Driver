%{
Run a 3D FEM simulation on a single subject
The driver will rotate the subject to have the origin at the bottom, posterior,
left corner of the torso. Positive y is towards the anterior side.
All units NEED to be in standard units to work.
Initial Version: 10/07/24 - Kyler Howard
Current Version: 01/30/26 - Kyler Howard
%}

close all
clearvars -except dataset_path old_dataset_path save_path old_save_path
clc
pause('on')

% Add paths in a way that works for macs as well
addpath(fullfile("Modified_OOEIT", "ForwardProblemSolvers"))
addpath(fullfile("Modified_OOEIT", "MiscClasses"))
addpath(fullfile("Utility_Functions"))

total_start_time = tic();

% ----------------------------------------------------------------------- %
%%                                Settings                                %
% ----------------------------------------------------------------------- %
% User settings
flags.do_pauses       = 0; % Decide to include pauses to check things or not
flags.solve_problem   = 1; % Decide if you want to setup (0), or fully solve (1)
flags.use_GE          = 1; % Decide if you want to use GE (1) or ACT5 (0) current patterns/conductivities
flags.do_parfor       = 1; % Decide if you want to paralize (1) or not (0)
flags.inject_current  = 1; % Decide if you want to inject ANY current (1) or only measure voltages (0)
flags.heart_BCs       = 0; % Decide if you want to include heart BCs (1) or not (0)
flags.save_heart_mesh = 0; % Decide if you want to generate and save a heart mesh (1) or not (0)
flags.do_beeps        = 0; % Decide if you want the code to beep after each simulation (1) or not (0)
flags.verbose         = 0; % Decide if you want to print status updates along the way (1) or not (0)

% Condition & permutation settings
    % Conditions are a cell array, where each cell contains its own condition
    % Condition is {max_inspiration, equal_vent, left_only, right_only, esoph_intubate, condition_name}
flags.conditions   = {{0.625, 1, 0, 0, 0, "Max_Insp"};...       Max Baby Inspiration
                      {0.375, 1, 0, 0, 0, "Min_Insp"};...       Max Baby Expiration
                      {0.500, 1, 0, 0, 0, "Mean_Insp"};...      Mean Inspiration
                    % {1.500, 1, 0, 0, 0, "Deep_Insp"},...      Deep Inspiration
                      {0.500, 0, 0, 1, 0, "Left_Intubate"};...  Left Bronchus Intubation
                      {0.500, 0, 1, 0, 0, "Right_Intubate"};... Right Bronchus Intubation
                      {0.000, 1, 0, 0, 1, "Esoph_Intubate"}}; % Esophageal Intubation
    % Permutations are a cell array, where each cell contains its own permutation settings
    % Perumutation is {num_perm, lung_range, esoph_range}
flags.permutations = {{2,  0.025, 0.000};... Max Baby Inspiration
                      {2,  0.025, 0.000};... Max Baby Expiration
                      {10, 0.125, 0.000};... Mean Inspiration
                    % {2,  0.100, 0.000};... Deep Inspiration
                      {5,  0.125, 0.000};... Left Bronchus Intubation
                      {5,  0.125, 0.000};... Right Bronchus Intubation
                      {5,  0.000, 0.125}}; % Esophageal Intubation

% Conductivity Settings
flags.set_complex       = 0; % Choice of complex (1) or real (0) conductivities
flags.const_body        = 0; % Decide if you want a solid/constant body (1) or not (0)
flags.do_conditions     = 1; % Decide if you want to run multiple conditions (1), or just the programed condition below (0)
    flags.max_inspiration   = 0.5; % Decide if the lungs should be at inspiration (1), expiration (0), or somewhere in-between
    flags.lung_range        = 0.25; % Decide what percentage of inspiration range you are okay with. Default is 0.25/25%
    flags.esoph_range       = 0.125; % Decide what percentage of inspiration range for the lungs for esoph intubation
    flags.equal_vent        = 1; % Decide if you want equal ventilation (1) or split (0)
    flags.left_only         = 0; % Decide if you want only ventilation on the left side (1) or not (0)
    flags.right_only        = 0; % Decide if you want only ventilation on the right side (1) or not (0)
    flags.esoph_intubate    = 0; % Decide if the esophagus is intubated (1) or not (0)
flags.permute_conds     = 1; % Decide if you want random conds (1) or not (0)

% Plot settings
flags.plot_slices     = 0; % Plot individual slices when going slice by slice
flags.plot_trachea    = 0; % Plotting of carina height & trachea orientation
flags.plot_electrodes = 0; % Plotting of electrode consturction
flags.plot_conds      = 0; % Plotting of conductivities
flags.plot_GTs        = 0; % Plot ground truth images
flags.plot_internal   = 0; % Plotting of internal nodes
flags.plot_volts      = 0; % Plotting of nodal voltages
flags.plot_heart      = 1; % Plot heart BCs
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

if isempty(gcp('nocreate')) && flags.solve_problem == 1 && flags.do_parfor == 1
    % Open a parallel pool
    parpool;
end

% ----------------------------------------------------------------------- %
%%                       Select The Dataset To Run                        %
% ----------------------------------------------------------------------- %

% Select the mesh you wish to run
if exist("dataset_path", "var") && ischar(dataset_path)
    % Use the previous mesh path
    old_dataset_path = dataset_path;
    [dataset_path]   = uigetdir(dataset_path, "Select Parent 'CTs' Folder With Data");
elseif exist("dataset_path", "var") && exist("old_dataset_path", "var")
    % User hit cancel. Use datset path from two attempts ago
    [dataset_path]   = uigetdir(old_dataset_path, "Select Parent 'CTs' Folder With Data");
else
    % First time running this
    [dataset_path]   = uigetdir(pwd, "Select Parent 'CTs' Folder With Data");
end

% Some user validation
if ischar(dataset_path) == 0
    error("Did not select a valid path")
end
if contains(dataset_path, "CT") == 0
    error("Did not select a path containing 'CTs'")
end

% ----------------------------------------------------------------------- %
%%                            Select Save Path                            %
% ----------------------------------------------------------------------- %

% Select the parent folder to save the file in
if exist("save_path", "var") && ischar(save_path)
    % Use the previous mesh path
    old_save_path = save_path;
    [save_path] = uigetdir(save_path, "Select Anatomical Atlas to Save Files In");
elseif exist("save_path", "var") && exist("old_save_path", "var")
    % User hit cancel. Use mesh path from two attempts ago
    [save_path] = uigetdir(old_save_path, "Select Anatomical Atlas to Save Files In");
else
    % First time running this
    [save_path] = uigetdir(dataset_path, "Select Anatomical Atlas to Save Files In");
end

% Some user validation
if ischar(save_path) == 0
    error("Did not select a valid path")
end
if contains(save_path, "Anatomical_Atlas") == 0
    error("Did not select a path containing 'Anatomical_Atlas'")
end

% ----------------------------------------------------------------------- %
%%                          Loop Through Dataset                          %
% ----------------------------------------------------------------------- %

dataset_contents = dir(dataset_path);
sbjs_solved = 0;
num_solves  = 0;
for sbj_i = 86:108%size(dataset_contents, 1)
    % Get the subject name, mesh path, and mesh name from the folder
    sbj_name = dataset_contents(sbj_i).name;
    msh_path = fullfile(dataset_path, sbj_name);
    msh_name = dir(fullfile(msh_path, "*Eroded.mat")).name;

    % Updating setting based on the mesh
    if contains(msh_name, 'NoBones')
        flags.are_bones = 0; % There are no bones
    else
        flags.are_bones = 1; % There are bones
    end

% ----------------------------------------------------------------------- %
%%                           Prep Saving Folders                          %
% ----------------------------------------------------------------------- %

    % Create the subject specific save path
    if contains(save_path, sbj_name) == 0
        sbj_save_path = fullfile(save_path, sbj_name);
        if isfolder(sbj_save_path) == 0
            mkdir(sbj_save_path)
        end
    end
    
    % Create the nested folders to save into
    electrode_save_path = fullfile(sbj_save_path, flags.E_type);
    if isfolder(electrode_save_path) == 0
        mkdir(electrode_save_path)
    end
    cond_save_path = fullfile(electrode_save_path, "conductivities");
    if isfolder(cond_save_path) == 0
        mkdir(cond_save_path)
    end
    volt_save_path = fullfile(electrode_save_path, "voltages");
    if isfolder(volt_save_path) == 0
        mkdir(volt_save_path)
    end
    % Loop through each condition
    for i = 1:numel(flags.conditions)
        condition_cond_save_path = fullfile(cond_save_path, flags.conditions{i}{6});
        if isfolder(condition_cond_save_path) == 0
            mkdir(condition_cond_save_path)
        end
        condition_volt_save_path = fullfile(volt_save_path, flags.conditions{i}{6});
        if isfolder(condition_volt_save_path) == 0
            mkdir(condition_volt_save_path)
        end
    end

% ----------------------------------------------------------------------- %
%%                         Running FEM3D Function                         %
% ----------------------------------------------------------------------- %
    
    sbj_start_time = tic;
    
    % RUN THE 3D FEM 
    fprintf("Running %s\n", sbj_name)
    FEM3D_Function(msh_path, msh_name, sbj_name, sbj_save_path, flags, noise)
    
    sbj_stop_time = toc(sbj_start_time);
    fprintf("\n   It took %.2f hours to solve the forward problem\n", sbj_stop_time / 3600)

    sbjs_solved = sbjs_solved + 1;
    num_solves  = num_solves + sum(cellfun(@(x) x{1}, flags.permutations));
end


total_stop_time = toc(total_start_time);
fprintf("\n\nIt took %.2f hours to solve %d subjects & %d unique solves\n", total_stop_time / 3600, sbjs_solved, num_solves)
if flags.do_beeps == 1
    beep
end
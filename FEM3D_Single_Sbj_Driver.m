%{
Run a 3D FEM simulation on a single subject
The driver will rotate the subject to have the origin at the bottom, posterior,
left corner of the torso. Positive y is towards the anterior side.
All units NEED to be in standard units to work.
Initial Version: 10/07/24 - Kyler Howard
Current Version: 01/30/26 - Kyler Howard
%}

all_fig = findall(0, "type", "figure");
close(all_fig)
clearvars -except msh_path old_msh_path save_path old_save_path
clc
pause("on")

% Add paths in a way that works for macs as well
addpath(fullfile("Modified_OOEIT", "ForwardProblemSolvers"))
addpath(fullfile("Modified_OOEIT", "MiscClasses"))
addpath(fullfile("Utility_Functions"))

% ----------------------------------------------------------------------- %
%%                                Settings                                %
% ----------------------------------------------------------------------- %
% User settings
flags.do_pauses       = 0; % Decide to include pauses to check things or not
flags.solve_problem   = 0; % Decide if you want to setup (0), or fully solve (1)
flags.use_GE          = 1; % Decide if you want to use GE (1) or ACT5 (0) current patterns/conductivities
flags.do_parfor       = 0; % Decide if you want to paralize (1) or not (0)
flags.inject_current  = 1; % Decide if you want to inject ANY current (1) or only measure voltages (0)
flags.heart_BCs       = 0; % Decide if you want to include heart BCs (1) or not (0)
flags.save_heart_mesh = 0; % Decide if you want to generate and save a heart mesh (1) or not (0)
flags.do_beeps        = 1; % Decide if you want the code to beep after each simulation (1) or not (0)
flags.verbose         = 1; % Decide if you want to print status updates along the way (1) or not (0)

% Video settings
flags.make_video  = 0;              % Decide if you want to make a video (1), or a single frame (0)
flags.breath_rate = 44;             % Breath rate in breaths per minute 44
flags.heart_rate  = 120;            % heart  rate in beats   per minute
flags.fps         = 28;             % Frame rate to reconstruct the video with
flags.insp_range  = [0.375, 0.625]; % Min and max inspiration percentages
% flags.insp_range  = [0, 1.5]; % Min and max inspiration percentages

% Condition & permutation settings
    % Conditions are a cell array, where each cell contains its own condition
    % Condition is {max_inspiration, equal_vent, left_only, right_only, esoph_intubate, condition_name}
flags.conditions   = {{0.500, 0, 0, 0, 0, "Reg_Intubate"};...   Regular Baby Inspiration
                    % {1.500, 1, 0, 0, 0, "Deep_Insp"},...      Deep Inspiration
                      {0.500, 0, 0, 1, 0, "Left_Intubate"};...  Left Bronchus Intubation
                      {0.500, 0, 1, 0, 0, "Right_Intubate"};... Right Bronchus Intubation
                      {0.000, 1, 0, 0, 1, "Esoph_Intubate"}}; % Esophageal Intubation
    % Permutations are a cell array, where each cell contains its own permutation settings
    % Perumutation is {num_perm, lung_range, esoph_range}
flags.permutations = {{10, 0.125, 0.000};... Regular Baby Inspiration
                    % {2,  0.100, 0.000};... Deep Inspiration
                      {5,  0.125, 0.000};... Left Bronchus Intubation
                      {5,  0.125, 0.000};... Right Bronchus Intubation
                      {5,  0.000, 0.125}}; % Esophageal Intubation

% Conductivity Settings
flags.set_complex       = 0; % Choice of complex (1) or real (0) conductivities
flags.const_body        = 0; % Decide if you want a solid/constant body (1) or not (0)
flags.do_conditions     = 0; % Decide if you want to run multiple conditions (1), or just the programed condition below (0)
    flags.max_inspiration   = 0.0; % Decide if the lungs should be at inspiration (1), expiration (0), or somewhere in-between. Also controls esoph in esoph_intubate.
    flags.cardiac_cycle     = 1;   % Decide if the heart should be at diastole (1), systole (0), or somewhere in-between
    flags.lung_range        = 0.25; % Decide what percentage of inspiration range you are okay with. Default is 0.25/25%
    flags.esoph_range       = 0.125; % Decide what percentage of inspiration range for the lungs for esoph intubation
    flags.equal_vent        = 1; % Decide if you want equal ventilation in each lung (1) or split (0)
    flags.left_only         = 0; % Decide if you want only ventilation on the left side (1) or not (0)
    flags.right_only        = 0; % Decide if you want only ventilation on the right side (1) or not (0)
    flags.esoph_intubate    = 0; % Decide if the esophagus is intubated (1) or not (0)
flags.permute_conds     = 0; % Decide if you want random conds (1) or not (0)

% Plot settings
flags.plot_slices     = 0; % Plot individual slices when going slice by slice
flags.plot_trachea    = 0; % Plotting of carina height & trachea orientation
flags.plot_electrodes = 1; % Plotting of electrode consturction
flags.plot_conds      = 0; % Plotting of conductivities
flags.plot_GTs        = 1; % Plot ground truth images
flags.plot_internal   = 0; % Plotting of internal nodes
flags.plot_volts      = 1; % Plotting of nodal voltages
flags.plot_heart      = 0; % Plot heart BCs
flags.plot_breath     = 1; % Plot the breathing and cardiac curves
flags.fixed_range     = 1; % Set GT plots to be a standard range


flags.CP_choice       = 1; % Choice of current pattern for patches
    % 1: Standard pattern
    % 2: 4x8 pattern
flags.E_choice        = 3; % Choice of Electrode configuration
    % 1: Large patch front back  (GE Patch)
    % 2: Small patch front back  (GE Patch)
    % 3: Two rows of large belts (GE Belt)
    % 4: Two rows of small belts (GE Belt)
    % 5: Custom electrodes

% Custom Electrode Settings
flags.E_type          = "belt";   % Choice between "patch" and "belt"
flags.E_shape         = "circle"; % Choice between "circle" and "rectangle"
flags.E_dia           = 30;       % Diameter of electrode in mm (for circle)
flags.E_width         = 22;       % Width  of electrode in mm (for rectangle)
flags.E_height        = 29;       % Height of electrode in mm (for rectangle)
flags.gap_width       = 46.675;   % Gap between electrodes horizontally in mm (edge-edge) (for patch) %2.5 / 46.675
flags.gap_height      = 32.875;   % Gap between electrodes vertically in mm (edge-edge) (for patch) %2.5 / 32.3875
flags.E_count         = [16];     % Number of electrodes per row (for belt), or matrix of how many rows and columns (for patch)
flags.equal_space     = 1;        % If the electrodes should be equally spaced (1) or start at the armpit and "rolled" on like GE (0)
flags.E_space         = 5;        % Edge-to-edge spacing between electrodes in mm for unequal belt spacing

% Save the correct electrode settings, not the custom ones when using the
% standard settings
[~, flags] = Construct_Electrode_Settings(flags);

noise = [0, 0, 0, 0];       % Noise and error parameters
    % noise_rel = err(1);
    % noise_abs = err(2);
    % e_systematic_rel = err(3);
    % e_systematic_abs = err(4);

if isempty(gcp("nocreate")) && flags.solve_problem == 1 && flags.do_parfor == 1
    % Open a parallel pool
    parpool;
end

% ----------------------------------------------------------------------- %
%%                         Select The Mesh To Run                         %
% ----------------------------------------------------------------------- %

% Select the mesh you wish to run
if exist("msh_path", "var") && ischar(msh_path)
    % Use the previous mesh path
    old_msh_path = msh_path;
    [msh_name, msh_path] = uigetfile(msh_path, "Select Mesh File");
elseif exist("msh_path", "var") && exist("old_msh_path", "var")
    % User hit cancel. Use mesh path from two attempts ago
    [msh_name, msh_path] = uigetfile(old_msh_path, "Select Mesh File");
else
    % First time running this
    [msh_name, msh_path] = uigetfile(pwd,"Select Mesh File");
end

% Some user validation
if ischar(msh_path) == 0
    error("Did not select a valid mesh path")
end
if contains(msh_name, ".mat") == 0
    error("Did not select a .mat file")
end
if contains(lower(msh_name), "mesh") == 0
    error("Did not select a valid mesh file")
end

% Extracting the subject name
parts    = split(msh_name, "_");
sbj_name = parts{1};

% Updating setting based on the mesh
if contains(msh_name, "NoBones")
    flags.are_bones = 0; % There are no bones
else
    flags.are_bones = 1; % There are bones
end

% ----------------------------------------------------------------------- %
%%                            Select Save Path                            %
% ----------------------------------------------------------------------- %

if flags.solve_problem == 1
    % Select the parent folder to save the file in
    if exist("save_path", "var") && ischar(save_path)
        % Use the previous mesh path
        old_save_path = save_path;
        [save_path] = uigetdir(save_path, "Select Parent Folder to Save Files In");
    elseif exist("save_path", "var") && exist("old_save_path", "var")
        % User hit cancel. Use mesh path from two attempts ago
        [save_path] = uigetdir(old_save_path, "Select Parent Folder to Save Files In");
    else
        % First time running this
        [save_path] = uigetdir(msh_path, "Select Parent Folder to Save Files In");
    end
    
    % Some user validation
    if ischar(save_path) == 0
        error("Did not select a valid save path")
    end
    
    % Check if the user isn't just saving to the results folder
    if contains(save_path, fullfile("OOEIT_FEM_Driver","Results")) ~= 1
        % Create the subject specific save path if the user didn't select it
        if contains(save_path, sbj_name) == 0
            sbj_save_path = fullfile(save_path, sbj_name);
            if isfolder(sbj_save_path) == 0
                mkdir(sbj_save_path)
            end
        else
            sbj_save_path = save_path;
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
    else
        sbj_save_path = save_path;
    end
else
    sbj_save_path = "";
end

% ----------------------------------------------------------------------- %
%%                         Running FEM3D Function                         %
% ----------------------------------------------------------------------- %

start_time = tic;

% RUN THE 3D FEM 
fprintf("Running %s\n", sbj_name)
n_bframes = FEM3D_Function(msh_path, msh_name, sbj_name, sbj_save_path, flags, noise);

stop_time = toc(start_time);

if flags.do_conditions == 1
    n_conditions = size(flags.conditions,1);
else
    n_conditions = 1;
end

num_solves  = n_conditions * n_bframes * ((sum(cellfun(@(x) x{1}, flags.permutations))-1)*flags.do_conditions*flags.permute_conds+1); % Adjusting for if we did permutations

fprintf("\n\nIt took %.2f hours to run %d unique solves\n", stop_time / 3600, num_solves)

if flags.do_beeps == 1
    beep
end

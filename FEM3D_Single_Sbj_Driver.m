%{
Run a 3D FEM simulation on a single subject
The driver will rotate the subject to have the origin at the bottom, posterior,
left corner of the torso. Positive y is towards the anterior side.
All units NEED to be in standard units to work.
Initial Version: 10/07/24 - Kyler Howard
Current Version: 01/30/26 - Kyler Howard
%}

close all
clearvars -except msh_path old_path
clc
pause('on')
start_time = tic;

% Add paths in a way that works for macs as well
addpath(fullfile("Modified_OOEIT", "ForwardProblemSolvers"))
addpath(fullfile("Modified_OOEIT", "MiscClasses"))
addpath(fullfile("Utility_Functions"))

% Select the mesh you wish to run
if exist("msh_path", "var") && ischar(msh_path)
    % Use the previous mesh path
    old_path = msh_path;
    [msh_name, msh_path] = uigetfile(msh_path, "Select Mesh File");
elseif exist("msh_path", "var") && exist("old_path", "var")
    % User hit cancel. Use mesh path from two attempts ago
    [msh_name, msh_path] = uigetfile(old_path, "Select Mesh File");
else
    % First time running this
    [msh_name, msh_path] = uigetfile("Select Mesh File");
end

% Some user validation
if ischar(msh_path) == 0
    error("Did not select a valid path")
end
if contains(msh_name, ".mat") == 0
    error("Did not select a .mat file")
end
if contains(lower(msh_name), "mesh") == 0
    error("Did not select a valid mesh file")
end

% Extracting the subject name
parts    = split(msh_name, '_');
sbj_name = parts{1};
fprintf("Running %s\n", sbj_name)

% ----------------------------------------------------------------------- %
%%                                Settings                                %
% ----------------------------------------------------------------------- %
% User settings
flags.do_pauses       = 0; % Decide to include pauses to check things or not
flags.solve_problem   = 1; % Decide if you want to setup (0), or fully solve (1)
flags.use_GE          = 1; % Decide if you want to use GE (1) or ACT5 (0) current patterns/conductivities
flags.do_parfor       = 0; % Decide if you want to paralize (1) or not (0)
flags.inject_current  = 1; % Decide if you want to inject ANY current (1) or only measure voltages (0)
flags.heart_BCs       = 0; % Decide if you want to include heart BCs (1) or not (0)
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
%%                         Running FEM3D Function                         %
% ----------------------------------------------------------------------- %
FEM3D_Function(msh_path, msh_name, sbj_name, flags, noise)

stop_time = toc(start_time);
fprintf("\n   It took %.2f hours to solve the forward problem\n", stop_time / 3600)

if flags.do_beeps == 1
    beep
end
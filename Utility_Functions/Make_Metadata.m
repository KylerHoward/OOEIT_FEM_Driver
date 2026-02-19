%{
Kyler Howard
Last Edit: 12/18/2025
Create metadata for each variable saved by the FEM3D Driver code
%}

function metadata = Make_Metadata(type)
    % Set file specific metadata
    if lower(type) == "volt"
        % Set variable descriptions
        metadata.descriptions.Umeas    = "Measured voltage at the center of each boundary condition. Saved as LxK (32x31).";
        metadata.descriptions.Uall     = "Measured voltage at every node in the mesh. Saved as nxK, number of nodes by currents.";
        metadata.descriptions.cur_pat  = "Applied current pattern, scaled if needed. Saved as LxK (32x31).";
        metadata.descriptions.perim_mm = "Perimeter of the body at the carina height (patch), or T5 height (belt).";
        metadata.descriptions.noise    = "Noise/Error vector. (1): relative noise, (2): absolute noise, (3): relative systematic error, (4): absolute systematic error.";

        % Set variable units
        metadata.units.Umeas    = "mV";
        metadata.units.Uall     = "mV";
        metadata.units.cur_pat  = "A";
        metadata.units.perim_mm = "mm";

    elseif lower(type) == "cond"
        % Set variable descriptions
        metadata.descriptions.sigma     = "Conductivity at every node in the mesh. Saved as nx1.";
        metadata.descriptions.sigma_GT  = "Cell array containing which sigma indices are in each ground truth image. Each cell is a binary nx1";
        metadata.descriptions.nodes     = "Matrix containing x/y/z coordinates for each node in the mesh. Saved as nx3";
        metadata.descriptions.E_connect = "Cell array containing each boundary condition and what surface faces are included. The values of the cell array are the indices of nodes. Each cell is mx3";

        % Set variable units
        metadata.units.sigma     = "S/m";
        metadata.units.sigma_GT  = "N/A, binary value";
        metadata.units.nodes     = "mm";
        metadata.units.E_connect = "N/A, node indices";
    end

    % Set flag descriptions
    metadata.descriptions.flags.do_pauses       = "Decide to include pauses to check things (1) or not (0). Default is 0.";
    metadata.descriptions.flags.solve_problem   = "Decide if you want to setup (0), or fully solve (1).";
    metadata.descriptions.flags.use_GE          = "Decide if you want to use GE (1) or ACT5 (0) current patterns/conductivities.";
    metadata.descriptions.flags.do_parfor       = "Decide if you want to paralize (1) or not (0) during solving.";
    metadata.descriptions.flags.inject_current  = "Decide if you want to inject ANY current (1) or only measure voltages (0). Default is 1.";
    metadata.descriptions.flags.heart_BCs       = "Decide if you want to include heart BCs (1) or not (0). Default is 0.";
    metadata.descriptions.flags.save_heart_mesh = "Decide if you want to generate and save a heart mesh (1) or not (0).";
    metadata.descriptions.flags.do_beeps        = "Decide if you want the code to beep after each simulation (1) or not (0).";
    metadata.descriptions.flags.verbose         = "Decide if you want to print status updates along the way (1) or not (0).";

    metadata.descriptions.flags.make_video      = "Decide if you want to make a video (1), or a single frame (0).";
    metadata.descriptions.flags.breath_rate     = "Breath rate in breaths per minute.";
    metadata.descriptions.flags.heart_rate      = "heart  rate in beats   per minute.";
    metadata.descriptions.flags.fps             = "Frame rate to reconstruct the video with.";
    metadata.descriptions.flags.insp_range      = "Min and max inspiration percentages.";

    metadata.descriptions.flags.conditions      = "Cell array, where each cell contains its own condition in the order of {max_inspiration, equal_vent, left_only, right_only, esoph_intubate, condition_name}.";
    metadata.descriptions.flags.permutations    = "Cell array, where each cell contains its own permutation settings in the order of {num_perm, lung_range, esoph_range}.";
    
    metadata.descriptions.flags.set_complex     = "Choice of complex (1) or real (0) conductivities. Default is 1.";
    metadata.descriptions.flags.const_body      = "Decide if you want a solid/constant body (1) or not (0). Default is 0.";
    metadata.descriptions.flags.do_conditions   = "Decide if you want to run multiple conditions (1), or just the programed condition below (0).";
    metadata.descriptions.flags.max_inspiration = "Decide if the lungs should be at inspiration (1), expiration (0), or somewhere in-between. Default is 0.5. Also controls esoph in esoph_intubate.";
    metadata.descriptions.flags.cardiac_cycle   = "Decide if the heart should be at diastole (1), systole (0), or somewhere in-between.";
    metadata.descriptions.flags.lung_range      = "Decide what percentage of inspiration range you are okay with. Default is 0.25/25.";
    metadata.descriptions.flags.esoph_range     = "Decide what percentage of inspiration range for the lungs for esoph intubation.";
    metadata.descriptions.flags.equal_vent      = "Decide if you want equal ventilation (1) or split (0).";
    metadata.descriptions.flags.left_only       = "Decide if you want only ventilation on the left side (1) or not (0).";
    metadata.descriptions.flags.right_only      = "Decide if you want only ventilation on the right side (1) or not (0).";
    metadata.descriptions.flags.esoph_intubate  = "Decide if the esophagus is intubated (1) or not (0).";
    metadata.descriptions.flags.permute_conds   = "Decide if you want random conds (1) or not (0).";

    metadata.descriptions.flags.plot_slices     = "Plot individual slices when going slice by slice.";
    metadata.descriptions.flags.plot_trachea    = "Plotting of carina height.";
    metadata.descriptions.flags.plot_electrodes = "Plotting of electrode consturction.";
    metadata.descriptions.flags.plot_conds      = "Plotting of conductivities.";
    metadata.descriptions.flags.plot_GTs        = "Plot ground truth images.";
    metadata.descriptions.flags.plot_internal   = "Plotting of internal nodes.";
    metadata.descriptions.flags.plot_volts      = "Plotting of nodal voltages.";
    metadata.descriptions.flags.plot_heart      = "Plot heart BCs.";
    metadata.descriptions.flags.plot_breath     = "Plot the breathing and cardiac curves.";
    metadata.descriptions.flags.fixed_range     = "Set GT plots to be a standard range.";

    metadata.descriptions.flags.CP_choice       = "Choice of current pattern for patches. 1 is the standard pattern, 2 is a 4x8 pattern for the patch.";
    metadata.descriptions.flags.E_choice        = "Choice of Electrode configuration. 1 is the large GE patch. 2 is the small GE patch. 3 is the large GE belt. 4 is the small GE belt. 5 is custom.";
    metadata.descriptions.flags.E_type          = "Choice between 'patch' and 'belt'.";
    metadata.descriptions.flags.E_shape         = "Choice between 'circle' and 'rectangle'.";
    metadata.descriptions.flags.E_dia           = "Diameter of electrode (for circle electrodes).";
    metadata.descriptions.flags.E_rad           = "Radius of electrodes (for circular electrodes).";
    metadata.descriptions.flags.E_width         = "Width  of electrode (for rectangle electrodes).";
    metadata.descriptions.flags.E_height        = "Height of electrode (for rectangle electrodes).";
    metadata.descriptions.flags.E_area          = "Area of electrodes.";
    metadata.descriptions.flags.gap_width       = "Gap between electrodes horizontally (edge-edge) (for patch electrodes).";
    metadata.descriptions.flags.gap_height      = "Gap between electrodes vertically (edge-edge) (for patch electrodes).";
    metadata.descriptions.flags.E_count         = "Number of electrodes per row (for belt), or matrix of how many rows and columns (for patch).";
    
    metadata.descriptions.flags.are_bones       = "Flag if the mesh has bones (1) or not (0).";
    metadata.descriptions.flags.zeta            = "Contact impedance on every boundary condition. Saved as Lx1.";
    metadata.descriptions.flags.breath_curve    = "Vector from 0:1 of max_inspiration values throughout a simulated 'video'.";
    metadata.descriptions.flags.heart_curve     = "Vector from 0:1 of cardiac_cycle values throughout a simulated 'video'.";

    % Set flag units
    metadata.units.flags.zeta        = "Ω m²";
    metadata.units.flags.E_dia       = "mm";
    metadata.units.flags.E_rad       = "mm";
    metadata.units.flags.E_width     = "mm";
    metadata.units.flags.E_height    = "mm";
    metadata.units.flags.gap_width   = "mm";
    metadata.units.flags.gap_height  = "mm";
    metadata.units.flags.E_area      = "mm²";
    metadata.units.flags.breath_rate = "breath/min";
    metadata.units.flags.heart_rate  = "beat/min";
    metadata.units.flags.fps         = "frame/s";

    % Store computer info
    metadata.computer_info.computer_type  = computer;
    metadata.computer_info.OS             = feature('GetOS');
    metadata.computer_info.n_cores        = feature('numCores');
    metadata.computer_info.n_processors   = str2double(getenv('NUMBER_OF_PROCESSORS')); 
    metadata.computer_info.CPU            = feature('GetCPU');
    metadata.computer_info.MATLAB_version = version;

    % Random metadata information
    metadata.timestamp = datetime('now', 'Format', 'm_dd_y HH:mm:ss');
    
    % Get the username depending on the operating system
    if ispc
        metadata.user = getenv('USERNAME');
    elseif isunix || ismac
        metadata.user = getenv('USER');
    end
end
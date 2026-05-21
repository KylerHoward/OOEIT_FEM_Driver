function [E, flags] = Construct_Electrode_Settings(flags)
    %{
    Construct the settings for GE electrode patch/belts or for a custom setup
    12/3/25 - Kyler Howard

    param: flags  - Various flags controlling plotting and other parameters

    return: E     - Structure containing all the settings of the electrodes
    return: flags - Updated various flags controlling plotting and other parameters
    %}

    % Constructing large GE patch
    L_square.type        = "patch";
    L_square.shape       = "rectangle";
    L_square.E_space     = NaN;              % mm
    L_square.E_width     = 10;               % mm
    L_square.E_height    = 10;               % mm
    L_square.E_dia       = NaN;
    L_square.E_rad       = NaN;
    L_square.E_area      = L_square.E_width * L_square.E_height;
    L_square.E_count     = [4,4];            % Electrodes per row and per column
    L_square.gap_width   = 2.5;              % mm (edge to edge)
    L_square.gap_height  = 2.5;              % mm (edge to edge)
    L_square.equal_space = NaN;
    
    % Constructing small GE patch
    S_square.type        = "patch";
    S_square.shape       = "rectangle";
    S_square.E_space     = NaN;              % mm
    S_square.E_width     = 7;                % mm
    S_square.E_height    = 7;                % mm
    S_square.E_dia       = NaN;
    S_square.E_rad       = NaN;
    S_square.E_area      = S_square.E_width * S_square.E_height;
    S_square.E_count     = [4,4];            % Electrodes per row and per column
    S_square.gap_width   = 2.5;              % mm (edge to edge)
    S_square.gap_height  = 2.5;              % mm (edge to edge)
    S_square.equal_space = NaN;
    
    % Constructing large GE belt
    L_belt.type        = "belt";
    L_belt.shape       = "circle";
    L_belt.E_space     = 6;                     % mm
    L_belt.E_dia       = 17;                    % mm
    L_belt.E_rad       = L_belt.E_dia / 2;      % mm
    L_belt.E_width     = NaN;
    L_belt.E_height    = NaN;
    L_belt.E_area      = pi * L_belt.E_rad^2;   % mm²
    L_belt.E_count     = 16;                    % Electrodes per row
    L_belt.gap_width   = NaN;
    L_belt.gap_height  = NaN;
    L_belt.equal_space = 1;
    
    % Constructing small GE belt
    S_belt.type        = "belt";
    S_belt.shape       = "circle";
    S_belt.E_space     = 5;                     % mm
    S_belt.E_dia       = 12;                    % mm
    S_belt.E_rad       = S_belt.E_dia / 2;      % mm
    S_belt.E_width     = NaN;
    S_belt.E_height    = NaN;
    S_belt.E_area      = pi * S_belt.E_rad^2;   % mm²
    S_belt.E_count     = 16;                    % Electrodes per row
    S_belt.gap_width   = NaN;
    S_belt.gap_height  = NaN;
    S_belt.equal_space = 1;

    % Constructing custom electrode setup
    E_custom.type  = flags.E_type;
    E_custom.shape = flags.E_shape;
    if E_custom.type == "patch"
        E_custom.E_count    = flags.E_count;
        E_custom.gap_width  = flags.gap_width;
        E_custom.gap_height = flags.gap_height;

        E_custom.E_space     = NaN;
        E_custom.equal_space = NaN;
    elseif E_custom.type == "belt"
        E_custom.E_count     = flags.E_count;
        E_custom.E_space     = flags.E_space;
        E_custom.equal_space = flags.equal_space;

        E_custom.gap_width  = NaN;
        E_custom.gap_height = NaN;
    end

    if E_custom.shape == "circle"
        E_custom.E_dia  = flags.E_dia;
        E_custom.E_rad  = flags.E_dia / 2;
        E_custom.E_area = pi * E_custom.E_rad^2;

        E_custom.E_width  = NaN;
        E_custom.E_height = NaN;
    elseif E_custom.shape == "rectangle"
        E_custom.E_width  = flags.E_width;
        E_custom.E_height = flags.E_height;
        E_custom.E_area   = E_custom.E_width * E_custom.E_height;

        E_custom.E_dia = NaN;
        E_custom.E_rad = NaN;
    end

    choices = {L_square, S_square, L_belt, S_belt, E_custom};
    E = choices{flags.E_choice};

    flags.E_type      = E.type;
    flags.E_shape     = E.shape;
    flags.E_width     = E.E_width;
    flags.E_height    = E.E_height;
    flags.E_dia       = E.E_dia;
    flags.E_rad       = E.E_rad;
    flags.E_area      = E.E_area;
    flags.E_count     = E.E_count;
    flags.gap_width   = E.gap_width;
    flags.gap_height  = E.gap_height;
    flags.E_space     = E.E_space;
    flags.equal_space = E.equal_space;
end
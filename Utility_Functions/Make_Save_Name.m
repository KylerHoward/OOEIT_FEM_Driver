function save_suffix = Make_Save_Name(condition_name, i_permutation, zeta, flags)

    save_suffix = "";
    if flags.do_conditions ~= 0
        save_suffix = sprintf("%s-%s%d", save_suffix, condition_name, i_permutation);
    end

    % Deal with the heart BC flags
    if flags.heart_BCs == 1
        save_suffix = sprintf("%s-HeartBCs%d", save_suffix);
    end
    if flags.inject_current == 0
        save_suffix = sprintf("%s-NoInjection", save_suffix);
    end

    % Add in the contact impedance
    save_suffix = sprintf("%s-z%g", save_suffix, zeta(1));

    % Check Constant Body
    if flags.const_body == 1
        save_suffix = sprintf("%s-SolidBody", save_suffix);
    elseif flags.do_conditions == 0
    % Otherwise, custom settings
        if flags.make_video == 1
            save_suffix = sprintf("%s-BreathCycle", save_suffix);
        else
            if flags.esoph_intubate == 1
                save_suffix = sprintf("%s-EsophInt", save_suffix);
            end
    
            save_suffix = sprintf("%s-%.2fInsp", save_suffix, flags.max_inspiration);
        end
        
        if flags.left_only == 1
            save_suffix = sprintf("%s-LeftOnly", save_suffix);
        elseif flags.right_only == 1
            save_suffix = sprintf("%s-RightOnly", save_suffix);
        end
    end
    
    % Electrode Settings
    if flags.E_choice == 1 
        save_suffix = sprintf("%s-LPatch", save_suffix);
    elseif flags.E_choice == 2
        save_suffix = sprintf("%s-SPatch", save_suffix);
    elseif flags.E_choice == 3
        save_suffix = sprintf("%s-LBelt", save_suffix);
    elseif flags.E_choice == 4
        save_suffix = sprintf("%s-SBelt", save_suffix);
    else
        % Custom Electrode Settings
        if flags.E_type == "belt"
            if flags.E_shape == "rectangle"
                save_suffix = sprintf("%s-Belt-RectE%dx%d", save_suffix, flags.E_width, flags.E_height);
            elseif flags.E_shape == "circle"
                save_suffix = sprintf("%s-Belt-CircleE%d", save_suffix, flags.E_dia);
            end
        elseif flags.E_type == "patch"
            if flags.E_shape == "rectangle"
                save_suffix = sprintf("%s-Patch-RectE%dx%d", save_suffix, flags.E_width, flags.E_height);
            elseif flags.E_shape == "circle"
                save_suffix = sprintf("%s-Patch-CircleE%d", save_suffix, flags.E_dia);
            end
        end
    end

    if flags.CP_choice == 2
        save_suffix = sprintf("%s-4x8", save_suffix);
    end
end
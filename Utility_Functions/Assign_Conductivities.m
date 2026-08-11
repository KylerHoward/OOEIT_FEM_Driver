    function sigma = Assign_Conductivities(nodes, connectivity, labels, lung_nodes, flags)
    %{
    Assign complex conductivities to each node
    8/20/24 - Kyler Howard

    param: nodes        - All nodes in the entire mesh (nx3)
    param: connectivity - Connectivity matrix for the entire mesh (mx4)
    param: labels       - Labels of each element for the entire mesh (mx1)
    param: flags        - Various flags controlling plotting and other parameters

    return: Sigma - Assigned complex conductivity of each node for the entire mesh at each frame (nNodes x nFrames)
    %}
    tic
%% ------------------------------- Setup -------------------------------- %

    % Set up Ventilation Lung Complex Conductivity settings. Averaging both frequencies from TFC
    nframes = size(flags.breath_curve, 2);
   
    % Set permitivity of free space
    eps0 = 8.8541878188e-12; % F/m

    % Set up system frequencies
    if flags.use_GE == 1
        freq = 10000; % Hz

        % Heart equations
        hc_full  = 0.72;
        hc_empty = 0.68;
        hc_m     = (hc_full - hc_empty) / (1-0);
        hc_eq    = hc_m * flags.heart_curve + hc_empty; % Updated to range from 0.55 - 0.75
        hc_range = 0.05 * ones([1, nframes]);

        hs_full  = 7.34e4; % 0.85 µF/m. Now 0.65
        hs_empty = 6.21e4; % 0.75 µF/m. Now 0.55
        hs_m     = (hs_full - hs_empty) / (1-0);
        hs_eq    = flags.heart_curve*hs_m + hs_empty;
        hs_range = 0.03 * ones([1, nframes]);

        % Lung equations
        lc_full  = 0.0932;
        lc_empty = 0.243;
        lc_m     = (lc_full - lc_empty) / (1-0);
        lc_eq    = lc_m * flags.breath_curve + lc_empty; % Linear range between empty/full
        lc_range = abs(lc_m * flags.lung_range) * ones([1, nframes]);

        ls_full  = 2.03e4; % 0.22 µF/m. Now 0.18
        ls_empty = 4.52e4; % 0.44 µF/m. Now 0.40
        ls_m     = (ls_full - ls_empty) / (1-0);
        ls_eq    = ls_m * flags.breath_curve + ls_empty; % Linear range between empty/full
        ls_range = (ls_m * flags.lung_range) * ones([1, nframes]); 
    
        % Dictionary of conductivities (S/m) from TFC. 
        % First value is the mean, second value is the allowable ± range
        conds.background  = [0;        0]      * ones([1, nframes]);   % 0.0,   0.0
        conds.soft_tissue = [0.400;    0.025]  * ones([1, nframes]);   % 0.4,   0.023
        conds.trachea     = [0.311;    0.0086] * ones([1, nframes]);   % 0.15,  0.0
        conds.bone        = [0.0204;   0.0003] * ones([1, nframes]);   % 0.05,  0.02
        conds.esophagus   = [0.530;    0.054]  * ones([1, nframes]);   % 0.164, 0.054
        conds.heart       = [hc_eq;    hc_range];                      % 0.66,  0.1
        conds.lung        = [lc_eq;    lc_range];
        conds.lung_tissue = [lc_empty; 0.002] * ones([1, nframes]);

        % Dictionary of complex conductivities (S/m). TFC gives multiples of Permittivity of free space
        % First value is the mean, second is the allowable ± range
        suscs.background  = 2*pi*freq*eps0*[0;        0]      * ones([1, nframes]);   % 0.0,  0.0
        suscs.soft_tissue = 2*pi*freq*eps0*[3.39e4;   5.65e3] * ones([1, nframes]);   % 0.2,  0.0
        suscs.trachea     = 2*pi*freq*eps0*[2.55e4;   2.26e3] * ones([1, nframes]);   % 0.05, 0.05
        suscs.bone        = 2*pi*freq*eps0*[5.20e2;   3.00e1] * ones([1, nframes]);   % 0.05, 0.05
        suscs.esophagus   = 2*pi*freq*eps0*[8.70e3;   3.00e2] * ones([1, nframes]);   % 0,    0
        suscs.heart       = 2*pi*freq*eps0*[hs_eq;    hs_range];                      % 0.4,  0.2
        suscs.lung        = 2*pi*freq*eps0*[ls_eq;    ls_range];                      % 0.4,  0.4
        suscs.lung_tissue = 2*pi*freq*eps0*[ls_empty; 1.15e4] * ones([1, nframes]);

    else % ACT5 System    
        freq = 93750; % Hz

        % Heart equations
        hc_full  = 0.72; % 0.75 or 0.72
        hc_empty = 0.68; % 0.55 or 0.68
        hc_m     = (hc_full - hc_empty) / (1-0);
        hc_eq    = flags.heart_curve*hc_m + hc_empty; % Updated to range from 0.55 - 0.75
        hc_range = 0.05 * ones([1, nframes]);

        hs_full  = 7.34e4; % 0.85 µF/m. Now 0.65
        hs_empty = 6.21e4; % 0.75 µF/m. Now 0.55
        hs_m     = (hs_full - hs_empty) / (1-0);
        hs_eq    = flags.heart_curve*hs_m + hs_empty;
        hs_range = 0.03 * ones([1, nframes]);

        % Lung equations
        lc_full  = 0.0932;
        lc_empty = 0.243;
        lc_m     = (lc_full - lc_empty) / (1-0);
        lc_eq    = flags.breath_curve*lc_m + lc_empty; % Linear range between empty/full
        lc_range = (lc_m * flags.lung_range) * ones([1, nframes]);

        ls_full  = 2.03e4; % 0.22 µF/m. Now 0.18
        ls_empty = 4.52e4; % 0.44 µF/m. Now 0.40
        ls_m     = (ls_full - ls_empty) / (1-0);
        ls_eq    = ls_m * flags.breath_curve + ls_empty; % Linear range between empty/full
        ls_range = (ls_m * flags.lung_range) * ones([1, nframes]); 
    
        % Dictionary of conductivity values from TFC. First value is the
        % mean, second value is the allowable ± range
        conds.background  = [0;        0]      * ones([1, nframes]);    % 0.0,   0.0
        conds.soft_tissue = [0.400;    0.025]  * ones([1, nframes]);    % 0.4,   0.023
        conds.trachea     = [0.337;    0.0086] * ones([1, nframes]);    % 0.15,  0.0
        conds.bone        = [0.0208;   0.0003] * ones([1, nframes]);    % 0.05,  0.02
        conds.esophagus   = [0.536;    0.054]  * ones([1, nframes]);    % 0.164, 0.054
        conds.heart       = [hc_eq;    hc_range];                       % 0.66 & 0.1
        conds.lung        = [lc_eq;    lc_range];
        conds.lung_tissue = [lc_empty; 0.002]  * ones([1, nframes]);

        % Dictionary of complex conductivities (S/m). TFC gives multiples of Permittivity of free space
        % First value is the mean, second is the allowable ± range
        suscs.background  = 2*pi*freq*eps0*[0;        0]      * ones([1, nframes]);   % 0.0,  0.0
        suscs.soft_tissue = 2*pi*freq*eps0*[3.39e4;   5.65e3] * ones([1, nframes]);   % 0.2,  0.0
        suscs.trachea     = 2*pi*freq*eps0*[2.55e4;   2.26e3] * ones([1, nframes]);   % 0.05, 0.05
        suscs.bone        = 2*pi*freq*eps0*[5.20e2;   3.00e1] * ones([1, nframes]);   % 0.05, 0.05
        suscs.esophagus   = 2*pi*freq*eps0*[8.70e3;   3.00e2] * ones([1, nframes]);   % 0,    0
        suscs.heart       = 2*pi*freq*eps0*[hs_eq;    hs_range];                      % 0.4,  0.2
        suscs.lung        = 2*pi*freq*eps0*[ls_eq;    ls_range];                      % 0.4,  0.4
        suscs.lung_tissue = 2*pi*freq*eps0*[ls_empty; 1.15e4] * ones([1, nframes]);
    end

    % Real Conductivities (S/m)
    background_cond  = conds.background(1,:);    background_cond_range  = conds.background(2,:);
    soft_tissue_cond = conds.soft_tissue(1,:);   soft_tissue_cond_range = conds.soft_tissue(2,:);
    trachea_cond     = conds.trachea(1,:);       trachea_cond_range     = conds.trachea(2,:);
    bone_cond        = conds.bone(1,:);          bone_cond_range        = conds.bone(2,:);
    heart_cond       = conds.heart(1,:);         heart_cond_range       = conds.heart(2,:);
    lung_cond        = conds.lung(1,:);          lung_cond_range        = conds.lung(2,:);
    lung_tissue_cond = conds.lung_tissue(1,:);   lung_tissue_cond_range = conds.lung_tissue(2,:);

    % Complex Conductivities (S/m)
    background_susc  = suscs.background(1,:);    background_susc_range  = suscs.background(2,:);
    soft_tissue_susc = suscs.soft_tissue(1,:);   soft_tissue_susc_range = suscs.soft_tissue(2,:);
    trachea_susc     = suscs.trachea(1,:);       trachea_susc_range     = suscs.trachea(2,:);
    bone_susc        = suscs.bone(1,:);          bone_susc_range        = suscs.bone(2,:);
    heart_susc       = suscs.heart(1,:);         heart_susc_range       = suscs.heart(2,:);
    lung_susc        = suscs.lung(1,:);          lung_susc_range        = suscs.lung(2,:);
    lung_tissue_susc = suscs.lung_tissue(1,:);   lung_tissue_susc_range = suscs.lung_tissue(2,:);

    % Real & complex values for if we are doing esophageal intubation
    if flags.esoph_intubate == 1
        esophagus_cond       = conds.lung(1,:);          esophagus_cond_range = lc_m * flags.esoph_range;
        esophagus_susc       = suscs.lung(1,:);          esophagus_susc_range = ls_m * flags.esoph_range;
        lung_cond            = conds.lung_tissue(1,:);   lung_cond_range      = conds.lung_tissue(2,:);
        lung_susc            = suscs.lung_tissue(1,:);   lung_susc_range      = suscs.lung_tissue(2,:);
    else
        esophagus_cond       = conds.esophagus(1,:);     esophagus_cond_range = conds.esophagus(2,:);
        esophagus_susc       = suscs.esophagus(1,:);     esophagus_susc_range = suscs.esophagus(2,:);
    end
    
    if flags.permute_conds == 1
        % Generate a number in [-1, 1] for each tissue type
        background_cond_pm  = 2*rand(1) - 1;
        background_susc_pm  = sign(background_cond_pm)*rand(1);
        left_lung_cond_pm   = 2*rand(1) - 1; 
        left_lung_susc_pm   = sign(left_lung_cond_pm)*rand(1);
        right_lung_cond_pm  = 2*rand(1) - 1; 
        right_lung_susc_pm  = sign(right_lung_cond_pm)*rand(1);
        trachea_cond_pm     = 2*rand(1) - 1; 
        trachea_susc_pm     = sign(trachea_cond_pm)*rand(1);
        soft_tissue_cond_pm = 2*rand(1) - 1; 
        soft_tissue_susc_pm = sign(soft_tissue_cond_pm)*rand(1);
        bone_cond_pm        = 2*rand(1) - 1; 
        bone_susc_pm        = sign(bone_cond_pm)*rand(1);
        esophagus_cond_pm   = 2*rand(1) - 1;
        esophagus_susc_pm   = sign(esophagus_cond_pm)*rand(1);
        heart_cond_pm       = -sign(left_lung_cond_pm)*rand(1);
        heart_susc_pm       = sign(heart_cond_pm)*rand(1);
    
        % Set complex conductivity values
        background_val  = background_cond + background_cond_range*background_cond_pm + ...
                          1i*(background_susc + background_susc_range*background_susc_pm);
        left_lung_val   = lung_cond + left_lung_cond_pm*lung_cond_range + ...
                          1i*(lung_susc + left_lung_susc_pm*lung_susc_range);
        if flags.equal_vent == 1
            right_lung_val = left_lung_val;
        elseif flags.equal_vent == 0
            right_lung_val  = lung_cond + right_lung_cond_pm*lung_cond_range + ...
                              1i*(lung_susc + right_lung_susc_pm*lung_susc_range);
        end
        trachea_val     = trachea_cond + trachea_cond_pm*trachea_cond_range + ...
                          1i*(trachea_susc + trachea_susc_pm*trachea_susc_range);
        soft_tissue_val = soft_tissue_cond + soft_tissue_cond_range*soft_tissue_cond_pm + ...
                          1i*(soft_tissue_susc + soft_tissue_susc_range*soft_tissue_susc_pm);
        if flags.left_only == 1
            right_lung_val = lung_tissue_cond;
            if flags.set_complex == 1
                right_lung_val = lung_tissue_cond + right_lung_cond_pm*lung_tissue_cond_range + ...
                                 1i*(lung_tissue_susc + right_lung_susc_pm*lung_tissue_susc_range);
            end
        elseif flags.right_only == 1
            left_lung_val = lung_tissue_cond;
            if flags.set_complex == 1
                 left_lung_val = lung_tissue_cond + left_lung_cond_pm*lung_tissue_cond_range + ...
                                 1i*(lung_tissue_susc + left_lung_susc_pm*lung_tissue_susc_range);
            end
        end
        bone_val        = bone_cond + bone_cond_pm*bone_cond_range + ...
                          1i*(bone_susc + bone_susc_pm*bone_susc_range);
        esophagus_val   = esophagus_cond + esophagus_cond_pm*esophagus_cond_range + ...
                          1i*(esophagus_susc + esophagus_susc_pm*esophagus_susc_range);
        heart_val       = heart_cond + heart_cond_range*heart_cond_pm + ...
                          1i*(heart_susc + heart_susc_range*heart_susc_pm);

        % Ensure each real value is positive
        left_lung_val  = make_real_positive(left_lung_val);
        right_lung_val = make_real_positive(right_lung_val);
    else
        % Set complex conductivity values
        background_val  = background_cond + 1i*background_susc;
        left_lung_val   = lung_cond + 1i*lung_susc;
        right_lung_val  = lung_cond + 1i*lung_susc;
        trachea_val     = trachea_cond + 1i*trachea_susc;
        soft_tissue_val = soft_tissue_cond + 1i*soft_tissue_susc;
        if flags.left_only == 1
            right_lung_val = lung_tissue_cond;
            if flags.set_complex == 1
                right_lung_val = lung_tissue_cond + 1i*lung_tissue_susc;
            end
        elseif flags.right_only == 1
            left_lung_val = lung_tissue_cond;
            if flags.set_complex == 1
                left_lung_val = lung_tissue_cond + 1i*lung_tissue_susc;
            end
        end
        bone_val        = bone_cond + 1i*bone_susc;
        esophagus_val   = esophagus_cond + 1i*esophagus_susc;
        heart_val       = heart_cond + 1i*heart_susc;
    end

    % Take just the real part if not doing complex
    if flags.set_complex == 0
        background_val  = real(background_val);
        left_lung_val   = real(left_lung_val);
        right_lung_val  = real(right_lung_val);
        trachea_val     = real(trachea_val);
        soft_tissue_val = real(soft_tissue_val);
        bone_val        = real(bone_val);
        esophagus_val   = real(esophagus_val);
        heart_val       = real(heart_val);
    end
    clear background_cond background_cond_range background_susc background_susc_range background_cond_pm background_susc_pm...
          lung_cond lung_cond_range lung_susc lung_susc_range left_lung_cond_pm left_lung_susc_pm right_lung_cond_pm right_lung_susc_pm...
          trachea_cond trachea_cond_range trachea_susc trachea_susc_range trachea_cond_pm trachea_susc_pm...
          soft_tissue_cond soft_tissue_cond_range soft_tissue_susc soft_tissue_susc_range soft_tissue_cond_pm soft_tissue_susc_pm...
          bone_cond bone_cond_range bone_susc bone_susc_range bone_cond_pm bone_susc_pm...
          esophagus_cond esophagus_cond_range esophagus_susc esophagus_susc_range esophagus_cond_pm esophagus_susc_pm...
          heart_cond heart_cond_range heart_susc heart_susc_range heart_cond_pm heart_susc_pm

    % Make a vector of conductivities in the same order as the labels
    if flags.are_bones == 1
        cond_vals = {background_val; {left_lung_val, right_lung_val}; trachea_val; soft_tissue_val; bone_val;...
                     esophagus_val; heart_val; background_val};
    else
        cond_vals = {background_val; {left_lung_val, right_lung_val}; trachea_val; soft_tissue_val;...
                     esophagus_val; heart_val; background_val};
    end
    % REMEMBER: LABELS ARE [0,8], NOT [1,9]. NEED A COND_VALS(LABEL + 1);

%% ----------------------------- Lung Prep ------------------------------ %
    split_correct = false;
    while split_correct == false
        [clusters, centroids] = kmeans(lung_nodes(:, [1,2]), 2, "Distance","cityblock");
        
        % Left and right are flipped
        if centroids(1,1) > centroids(2,1)
            clusters(clusters==1) = 3;
            clusters(clusters==2) = 1;
            clusters(clusters==3) = 2;

            centroids(3,:) = centroids(1,:);
            centroids(1,:) = centroids(2,:);
            centroids(2,:) = centroids(3,:);
            centroids(3,:) = [];
        end
    
        x_diff = abs(centroids(1,1) - centroids(2,1));
        y_diff = abs(centroids(1,2) - centroids(2,2));
        if x_diff > y_diff
            if flags.verbose == 1
                fprintf("      Split lungs correctly\n")
            end
            split_correct = true;
        else
            if flags.verbose == 1
                fprintf("      Split lungs failed\n")
            end
        end
    end
        
    if flags.plot_conds == 1
        figure()
            hold on
            scatter3(lung_nodes(clusters==1,1),lung_nodes(clusters==1,2),lung_nodes(clusters==1,3),"b")
            scatter3(lung_nodes(clusters==2,1),lung_nodes(clusters==2,2),lung_nodes(clusters==2,3),"red")
            xlabel("x")
            ylabel("y")
            zlabel("z")
            legend("Left (1)", "Right (2)", "location", "southoutside")
            title("xy projected - cityblock distance")
    end


%% ----------------------------- Assigning ------------------------------ %

    % determine counts for each node
    [~,~,idx] = unique(connectivity);
    counts    = accumarray(idx(:), 1);

    % Initialize the conductivity matrix to be NaNs and be large enough for
    % multiple elements to use the same node
    sigma = cell(size(nodes,1), nframes);
    for node = 1:size(nodes,1) % Go to the maximum node used
        [sigma{node,:}] = deal(nan([1, counts(node)], "like", cond_vals{3}));
    end

    % Cycle through elements
    for element = 1:size(connectivity,1)
        label = labels(element);
        % Cycle through the four nodes that make up each element
        for i = 1:4
            % Extract the node number
            node = connectivity(element, i);
            
            % Find the first column that is not a NaN
            col = 1;
            % while ~isnan(sigma(node, col))
            while ~isnan(sigma{node,1}(col))
                col = col + 1;
            end
            
            % Assign the conductivity value for constant tissues
            if label  ~= 1
                for c = 1:nframes
                    sigma{node,c}(col) = cond_vals{label + 1}(c);
                end

            % Assign the conductivity to each lung
            else
                % Find which cluster the node is in (1 = left, 2 = right)
                dist       = sum(abs(nodes(node,[1,2]) - centroids(:,:)), 2);
                [~, index] = min(dist);
                for c = 1:nframes
                    sigma{node,c}(col) = cond_vals{label + 1}{index}(c);
                end
            end

        end        
    end

    % Take a row wise average of the conductivities
    for r = 1:size(sigma,1)
        for c = 1:nframes
            sigma{r,c} = mode(sigma{r,c}, 2);
        end
    end
    sigma = cell2mat(sigma);

    % Set constant if wanted
    if flags.const_body == 1
        % Find all nodes that are not background
        sigma(~(real(sigma) == 0 & imag(sigma) == 0)) = cond_vals{4};
    end

%% ----------------------------- Plotting ------------------------------- %
    if flags.plot_conds == 1
        figure()
        if flags.set_complex == 1
            ax1 = subplot(1,2,1);
                scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 10, real(sigma(:,1)))
                title("Conductivity")
                xlabel("x (mm)")
                ylabel("y (mm)")
                zlabel("z (mm)")
                c = colorbar;
                c.Label.String = "S/m";
                c.Location     = "southoutside";
    
            ax2 = subplot(1,2,2);
                scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 10, imag(sigma(:,1)))
                title("Susceptivity")
                xlabel("x (mm)")
                ylabel("y (mm)")
                zlabel("z (mm)")
                c = colorbar;
                c.Label.String = "S/m";
                c.Location     = "southoutside";

            % Add a slider for x-axis min limit
            uicontrol("Style", "slider", "Min", min(nodes(:,1)), "Max", max(nodes(:,1))-1, "Value", 0, ...
                      "Position", [210, 30, 150, 20], ...
                      "Callback", @(src, event) update_xlim([ax1, ax2], src.Value, max(nodes(:,1))));
    
            % Add a slider for y-axis min limit
            uicontrol("Style", "slider", "Min", min(nodes(:,2)), "Max", max(nodes(:,2))-1, "Value", 0, ...
                      "Position", [210, 10, 150, 20], ...
                      "Callback", @(src, event) update_ylim([ax1, ax2], src.Value, max(nodes(:,1))));
        else
            scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 10, real(sigma(:,1)))
                title("Conductivity")
                xlabel("x (mm)")
                ylabel("y (mm)")
                zlabel("z (mm)")
                c = colorbar;
                c.Label.String = "S/m";
                c.Location     = "southoutside";

            % Add a slider for x-axis min limit
            uicontrol("Style", "slider", "Min", min(nodes(:,1)), "Max", max(nodes(:,1))-1, "Value", 0, ...
                      "Position", [210, 30, 150, 20], ...
                      "Callback", @(src, event) update_xlim(gca, src.Value, max(nodes(:,1))));
    
            % Add a slider for y-axis min limit
            uicontrol("Style", "slider", "Min", min(nodes(:,2)), "Max", max(nodes(:,2))-1, "Value", 0, ...
                      "Position", [210, 10, 150, 20], ...
                      "Callback", @(src, event) update_ylim(gca, src.Value, max(nodes(:,1))));
        end

        if flags.do_pauses == 1
            fprintf("      Press any key if correct\n")
            pause
        end
    end
end


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

function value = make_real_positive(value)
    % Ensure that only the real part gets set to zero
    if real(value) < 0
        value = 0 + imag(value);
    end
end
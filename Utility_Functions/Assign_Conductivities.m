function sigma = Assign_Conductivities(nodes, connectivity, labels, lung_nodes, flags)
    %{
    Assign complex conductivities to each node
    8/20/24 - Kyler Howard

    param: nodes        - All nodes in the entire mesh (nx3)
    param: connectivity - Connectivity matrix for the entire mesh (mx4)
    param: labels       - Labels of each element for the entire mesh (mx1)
    param: flags        - Various flags controlling plotting and other parameters

    return: Sigma - Assigned complex conductivity of each node for the entire mesh (nx1)
    %}
    tic
%% ------------------------------- Setup -------------------------------- %

    % S/m
    % Set up Ventilation Lung Conductivity settings
    if flags.use_GE == 1 % Freq = 10 kHz
        background_cond  = 0.0;    background_cond_range  = 0.0;     % 0.0,   0.0
        trachea_cond     = 0.311;  trachea_cond_range     = 0.0086;  % 0.15,  0.0
        soft_tissue_cond = 0.4;    soft_tissue_cond_range = 0.023;   % 0.3,   0.0       % 0.4, 0.023
        bone_cond        = 0.0204; bone_cond_range        = 0.0003;  % 0.05,  0.02
        % heart_cond       = 0.66;   heart_cond_range       = 0.1;     % 0.4,   0.1
        % KH: Updated 2/18/26 based on TFC and the cardiac cycle
        % heart_cond       = (50+70*flags.cardiac_cycle)/120*0.6 + (70*flags.cardiac_cycle)/120*0.153;
        heart_cond       = flags.cardiac_cycle*.2 + .55; % Updated to range from 0.55 - 0.75
        heart_cond_range = 0.05;
        % lung_cond        = 0.168;   lung_cond_range        = 0.075;    % 0.175, 0.125
        % KH: Updated 2/13/25 based on TFC
        lung_m           = -0.1498;
        lung_b           =  0.243;
        lung_cond        = lung_m * flags.max_inspiration + lung_b;  % Linear range between max/min
        lung_cond_range  = abs((lung_m * 0 + lung_b) - (lung_m * flags.lung_range + lung_b)); 
        if flags.esoph_intubate == 1
            esophagus_cond   = lung_m * flags.max_inspiration + lung_b;  
            esophagus_cond_range   = abs((lung_m * 0 + lung_b) - (lung_m * flags.esoph_range + lung_b));   % 0.164, 0.054
            lung_cond        = 0.243;  lung_cond_range = 0;
        else
            esophagus_cond   = 0.530;  esophagus_cond_range   = 0.054;   % 0.164, 0.054
        end
        lung_tissue_cond = 0.243;
    
        % Set up Complex Conductivity settings
        background_susc  = 0.0;    background_susc_range  = 0.0;     % 0.0,   0.0
        lung_susc        = 0.4;    lung_susc_range        = 0.4;     % 0.4,   0.4
        trachea_susc     = 0.05;   trachea_susc_range     = 0.05;    % 0.05,  0.05
        soft_tissue_susc = 0.2;    soft_tissue_susc_range = 0.0;     % 0.2,   0.0
        bone_susc        = 0.05;   bone_susc_range        = 0.05;    % 0.05,  0.05
        esophagus_susc   = 0.0;    esophagus_susc_range   = 0.0;     % 0.0,   0.0
        heart_susc       = 0.4;    heart_susc_range       = 0.2;     % 0.4,   0.2

    else % Freq = 93 kHz
        error("KYLER YOU SHOULD REALLY FIX THIS")
        if flags.verbose == 1
            fprintf("      All these values are wrong. Doing this for Chris\n")
        end
        background_cond  = 0.0;    background_cond_range  = 0.0;     % 0.0,   0.0
        trachea_cond     = 0.311;  trachea_cond_range     = 0.0086;  % 0.15,  0.0
        soft_tissue_cond = 0.4;    soft_tissue_cond_range = 0.023;   % 0.3,   0.0       % 0.2254, 0.023
        bone_cond        = 0.0204; bone_cond_range        = 0.0003;  % 0.05,  0.02
        heart_cond       = 0.66;   heart_cond_range       = 0.1;     % 0.4,   0.1
        % lung_cond        = 0.168;   lung_cond_range        = 0.075;    % 0.175, 0.125
        % KH: Updated 2/13/25 based on TFC
        lung_cond        = -0.1498 * flags.max_inspiration + 0.243;  % Linear range between max/min
        lung_cond_range  = 0.125;
        if flags.esoph_intubate == 1
            esophagus_cond   = 0.164;  esophagus_cond_range   = 0.054;   % 0.164, 0.054
        else
            esophagus_cond   = 0.530;  esophagus_cond_range   = 0.054;   % 0.164, 0.054
        end
    
        % Set up Complex Conductivity settings
        background_susc  = 0.0;    background_susc_range  = 0.0;     % 0.0,   0.0
        lung_susc        = 0.4;    lung_susc_range        = 0.4;     % 0.4,   0.4
        trachea_susc     = 0.05;   trachea_susc_range     = 0.05;    % 0.05,  0.05
        soft_tissue_susc = 0.2;    soft_tissue_susc_range = 0.0;     % 0.2,   0.0
        bone_susc        = 0.05;   bone_susc_range        = 0.05;    % 0.05,  0.05
        esophagus_susc   = 0.0;    esophagus_susc_range   = 0.0;     % 0.0,   0.0
        heart_susc       = 0.4;    heart_susc_range       = 0.2;     % 0.4,   0.2
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
                error("KYLER FIX THIS")
            end
        elseif flags.right_only == 1
            left_lung_val = lung_tissue_cond;
            if flags.set_complex == 1
                error("KYLER FIX THIS")
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
                error("KYLER FIX THIS")
            end
        elseif flags.right_only == 1
            left_lung_val = lung_tissue_cond;
            if flags.set_complex == 1
                error("KYLER FIX THIS")
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
        [clusters, centroids] = kmeans(lung_nodes(:, [1,2]), 2, 'Distance','cityblock');
        
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
            scatter3(lung_nodes(clusters==1,1),lung_nodes(clusters==1,2),lung_nodes(clusters==1,3),'b')
            scatter3(lung_nodes(clusters==2,1),lung_nodes(clusters==2,2),lung_nodes(clusters==2,3),'red')
            xlabel('x')
            ylabel('y')
            zlabel('z')
            legend("Left (1)", "Right (2)", 'location', 'southoutside')
            title('xy projected - cityblock distance')
    end


%% ----------------------------- Assigning ------------------------------ %

    % determine counts for each node
    [~,~,idx] = unique(connectivity);
    counts    = accumarray(idx(:), 1);

    % Initialize the conductivity matrix to be NaNs and be large enough for
    % multiple elements to use the same node
    sigma = cell(size(nodes,1), 1);
    for node = 1:size(nodes,1)
        sigma{node} = nan([1, counts(node)], 'like', cond_vals{3});
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
            while ~isnan(sigma{node}(col))
                col = col + 1;
            end
            
            % Assign the conductivity value for constant tissues
            if label  ~= 1
                sigma{node}(col) = cond_vals{label + 1};

            % Assign the conductivity to each lung
            else
                % Find which cluster the node is in (1 = left, 2 = right)
                dist       = sum(abs(nodes(node,[1,2]) - centroids(:,:)), 2);
                [~, index] = min(dist);
                sigma{node}(col) = cond_vals{label + 1}{index};
            end

        end        
    end

    % Take a row wise average of the conductivities
    sigma = cellfun(@(x) mode(x, 2), sigma);

    % Check for weird gaps from Cleaver
    % zero_nodes   = nodes(real(sigma) == 0,:);
    % zero_indices = Find_Internal_Nodes(zero_nodes, flags);
    % [~, zero_indices, ~] = intersect(nodes, zero_nodes(zero_indices,:),'rows');
    % sigma(zero_indices)  = soft_tissue_val;

    % Set constant if wanted
    if flags.const_body == 1
        % Find all nodes that are not background
        sigma(~(real(sigma) == 0 & imag(sigma) == 0)) = cond_vals{4};
    end

    cond_time = toc;
    if flags.verbose == 1
        fprintf("      It took %.2f seconds to assign values\n", cond_time)
    end

%% ----------------------------- Plotting ------------------------------- %
    if flags.plot_conds == 1
        figure()
        if flags.set_complex == 1
            ax1 = subplot(1,2,1);
                scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 10, real(sigma))
                title("Conductivity")
                xlabel("x (mm)")
                ylabel("y (mm)")
                zlabel("z (mm)")
                c = colorbar;
                c.Label.String = "S/m";
                c.Location     = 'southoutside';
    
            ax2 = subplot(1,2,2);
                scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 10, imag(sigma))
                title("Susceptivity")
                xlabel("x (mm)")
                ylabel("y (mm)")
                zlabel("z (mm)")
                c = colorbar;
                c.Label.String = "S/m";
                c.Location     = 'southoutside';

            % Add a slider for x-axis min limit
            uicontrol('Style', 'slider', 'Min', min(nodes(:,1)), 'Max', max(nodes(:,1))-1, 'Value', 0, ...
                      'Position', [210, 30, 150, 20], ...
                      'Callback', @(src, event) update_xlim([ax1, ax2], src.Value, max(nodes(:,1))));
    
            % Add a slider for y-axis min limit
            uicontrol('Style', 'slider', 'Min', min(nodes(:,2)), 'Max', max(nodes(:,2))-1, 'Value', 0, ...
                      'Position', [210, 10, 150, 20], ...
                      'Callback', @(src, event) update_ylim([ax1, ax2], src.Value, max(nodes(:,1))));
        else
            scatter3(nodes(:,1), nodes(:,2), nodes(:,3), 10, real(sigma))
                title("Conductivity")
                xlabel("x (mm)")
                ylabel("y (mm)")
                zlabel("z (mm)")
                c = colorbar;
                c.Label.String = "S/m";
                c.Location     = 'southoutside';

            % Add a slider for x-axis min limit
            uicontrol('Style', 'slider', 'Min', min(nodes(:,1)), 'Max', max(nodes(:,1))-1, 'Value', 0, ...
                      'Position', [210, 30, 150, 20], ...
                      'Callback', @(src, event) update_xlim(gca, src.Value, max(nodes(:,1))));
    
            % Add a slider for y-axis min limit
            uicontrol('Style', 'slider', 'Min', min(nodes(:,2)), 'Max', max(nodes(:,2))-1, 'Value', 0, ...
                      'Position', [210, 10, 150, 20], ...
                      'Callback', @(src, event) update_ylim(gca, src.Value, max(nodes(:,1))));
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
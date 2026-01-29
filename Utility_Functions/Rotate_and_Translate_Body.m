function [nodes, sbj_info] = Rotate_and_Translate_Body(nodes, organ_connects, carina_height, T5_height, T8_height, flags)
    %{
    Assess if a torso is rotated or shifted, and put it back in a standard
    orientation where:
        Anterior is positive y
        Superior is positive z
        Right    is positive x
    8/22/25 - Kyler Howard
    
    param: nodes          - Tetrahedron node list, n by 3, where n is the number of nodes
    param: organ_connects - Cell array of connectivity for each organ, n by 1, where n is the number of organs
                            Each cell is n by 4, where n is the number of elements
    param: carina_height  - Carina height from the excel file
    param: T5_height      - T5 height from the excel file
    param: T8_height      - T8 height from the excel file
    param: flags          - Various flags controlling plotting and other parameters
    
    return: nodes    - Rotated and translated nodes
    return: sbj_info - Structure containing translated carina, T5, and T8 heights
    %}

    % Create the organ numbering for reference 
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
        background  = 1;
        lung        = 2;
        trachea     = 3;
        soft_tissue = 4;
        esophagus   = 5;
        heart       = 6;
        external    = 7;
        % background = 1;
        % soft_tissue = 2;
        % external = 3;
    end
    
    % Determine if the body is upside down
    [body_nodes, ~]     = Get_Tet_Nodes(nodes, organ_connects{soft_tissue});
    [trachea_nodes, ~]  = Get_Tet_Nodes(nodes, organ_connects{trachea});
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
    end
    
    % Determine if the body is not facing towards y (anterior) side
    [heart_nodes, ~]   = Get_Tet_Nodes(nodes, organ_connects{heart});
    [body_nodes, ~]    = Get_Tet_Nodes(nodes, organ_connects{soft_tissue});
    [trachea_nodes, ~] = Get_Tet_Nodes(nodes, organ_connects{trachea});
    
    heart_center       = range(heart_nodes)/2 + min(heart_nodes);
    body_center        = range(body_nodes)/2 + min(body_nodes);
    trachea_center     = range(trachea_nodes)/2 + min(trachea_nodes);
    
    heart_diff = (heart_center   - body_center);
    trach_diff = (trachea_center - body_center);
    heart_range = range(heart_nodes);
    if ~(heart_diff(2) > 0.05*heart_range(2)) || ~(-trach_diff(2) > 0.05*heart_range(2)) % Check if the heart/trachea is on the wrong side
        % Facing -y (posterior) 
        if (-heart_diff(2) > 0.05*heart_range(2)) || (trach_diff(2) > 0.05*heart_range(2)) 
            theta = pi;
        % Facing x (right)
        elseif (heart_diff(1) > 0.05*heart_range(1)) || (-trach_diff(1) > 0.05*heart_range(1)) 
            theta = -pi / 2;
        % Facing -x (left)
        elseif (-heart_diff(1) > 0.05*heart_range(1)) || (trach_diff(1) > 0.05*heart_range(1)) 
            theta = pi / 2;
        end

        fprintf("Rotating the Chest Forward\n")
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

end

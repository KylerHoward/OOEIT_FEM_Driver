function [nodes, sbj_info] = Rotate_and_Translate_Body(nodes, labels, organ_connects, carina_height, T5_height, T8_height, flags)
    %{
    Assess if a torso is rotated or shifted, and put it back in a standard
    orientation where:
        Anterior is positive y
        Superior is positive z
        Right    is positive x
    8/22/25  - Kyler Howard
    12/11/25 - Edited for GMSH Kyler Howard
    
    param: nodes          - Tetrahedron node list, n by 3, where n is the number of nodes
    param: labels         - Tetrahedron label list, n by 1, where n is the number of elements
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
        lung        = 1;
        trachea     = 2;
        soft_tissue = 3;
        bone        = 4;
        if length(unique(labels)) == 6
            esophagus   = 5;
            heart       = 6;
        elseif length(unique(labels)) == 5
            heart = 5;
        end
    else
        lung        = 1;
        trachea     = 2;
        soft_tissue = 3;
        if length(unique(labels)) == 5
            esophagus   = 4;
            heart       = 5;
        elseif length(unique(labels)) == 4
            heart = 4;
        end
    end

% ----------------------------------------------------------------------- %
%%                            Sit Body Upright                            %
% ----------------------------------------------------------------------- %
    % Determine if the body is on its side
    [trachea_nodes, ~] = Get_Tet_Nodes(nodes, organ_connects{trachea});
    [~, min_index]     = min(range(trachea_nodes));
    [~, max_index]     = max(range(trachea_nodes));
    % Check if the minimm range of the trachea is on the z axis -> laying down
    % KH: 1/8/26. Swapping to a max index not being z instead
    if min_index == 3 || max_index ~= 3
        if flags.verbose == 1
            fprintf("   Sitting Body Upright\n")
        end
        theta = pi/2;
        rotationMatrix = [1, 0,           0;...
                          0, cos(theta), -sin(theta);...
                          0, sin(theta), cos(theta)];
    
        % Rotate and shift the body
        nodes      = nodes*rotationMatrix;
        nodes(:,2) = nodes(:,2) + abs(min(nodes(:,2)));
    end
    
% ----------------------------------------------------------------------- %
%%                           Flip Rightside Up                            %
% ----------------------------------------------------------------------- %
    % Determine if the body is upside down
    [body_nodes, ~]     = Get_Tet_Nodes(nodes, organ_connects{soft_tissue});
    [trachea_nodes, ~]  = Get_Tet_Nodes(nodes, organ_connects{trachea});
    % Check if the difference between the top of the trachea and the top of the body is less than 2 mm
    % Using the percentage of height doesn't work for short subjects
    if abs(max(trachea_nodes(:,3)) - max(body_nodes(:,3))) > 2
        if flags.verbose == 1
            fprintf("   Rotating Body Rightside Up\n")
        end
        theta = pi;
        rotationMatrix = [cos(theta), 0, -sin(theta);...
                          0,          1,  0;...
                          sin(theta), 0,  cos(theta)];
    
        % Rotate and shift the body
        nodes      = nodes*rotationMatrix;
        nodes(:,1) = nodes(:,1) * -1;
        nodes(:,3) = nodes(:,3) + abs(min(nodes(:,3)));
    end
    
% ----------------------------------------------------------------------- %
%%                              Face Forward                              %
% ----------------------------------------------------------------------- %
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

        if flags.verbose == 1
            fprintf("   Rotating the Chest Forward\n")
        end
        rotationMatrix = [cos(theta), -sin(theta), 0;...
                          sin(theta),  cos(theta), 0;...
                          0,           0,          1];
    
        % Rotate and shift the body
        nodes      = nodes*rotationMatrix;
        nodes(:,1) = nodes(:,1) + abs(min(nodes(:,1)));
        nodes(:,2) = nodes(:,2) + abs(min(nodes(:,2)));
    end
    
% % ----------------------------------------------------------------------- %
% %%                             Orient Trachea                             %
% % ----------------------------------------------------------------------- %       
%     % Rotate to make trachea straight up
%     if flags.verbose == 1
%         fprintf("   Lining the Trachea Vertical\n")
%     end
% 
%     % Check the current trachea angle
%     [trachea_nodes, ~] = Get_Tet_Nodes(nodes, organ_connects{trachea});
%     [body_nodes, ~]    = Get_Tet_Nodes(nodes, organ_connects{soft_tissue});
%     mean_body = mean(body_nodes);
%     trach_nodes = [trachea_nodes(trachea_nodes(:,3)>=carina_height,1), trachea_nodes(trachea_nodes(:,3)>=carina_height,3)];
%     mean_bot = mean(trach_nodes(trach_nodes(:,2)<min(trach_nodes(:,2))+0.1,:));
%     mean_top = mean(trach_nodes(trach_nodes(:,2)>max(trach_nodes(:,2))-0.1,:));
%     endpoints = [mean_bot; mean_top];
%     slope = (mean_bot(2) - mean_top(2)) / (mean_bot(1) - mean_top(1));
%     carina = [mean_bot, carina_height];
%     T5     = [mean_body(1), mean_body(2), T5_height];
%     T8     = [mean_body(1), mean_body(2), T8_height];
% 
%     if flags.plot_trachea == 1
%         figure; 
%             hold on
%             scatter(trachea_nodes(:,1), trachea_nodes(:,3))
%             scatter(mean_bot(1), mean_bot(2), "r", "filled")
%             scatter(mean_top(1), mean_top(2), "r", "filled")
%             plot(endpoints(:,1), endpoints(:,2),"r")
%             title(sprintf("Slope of %.2f\n", slope))
%     end
% 
%     if sign(slope) == 1
%         theta = -pi/50;
%     else
%         theta = pi/50;
%     end
%     rotationMatrix = [cos(theta), 0, -sin(theta);...
%                       0,          1,  0;...
%                       sin(theta), 0,  cos(theta)];
% 
%     % Rotate and shift the body
%     nodes  = nodes*rotationMatrix;
%     carina = carina*rotationMatrix;
%     T5     = T5*rotationMatrix;
%     T8     = T8*rotationMatrix;
% 
%     slope_new   = 0;
%     sign_change = 1;
%     while abs(slope_new) < 150
%         [trachea_nodes, ~] = Get_Tet_Nodes(nodes, organ_connects{trachea});
%         [body_nodes, ~]    = Get_Tet_Nodes(nodes, organ_connects{soft_tissue});
%         mean_body = mean(body_nodes);
%         trach_nodes = [trachea_nodes(trachea_nodes(:,3)>=carina(3),1), trachea_nodes(trachea_nodes(:,3)>=carina(3),3)];
%         mean_bot = mean(trach_nodes(trach_nodes(:,2)<min(trach_nodes(:,2))+0.1,:));
%         mean_top = mean(trach_nodes(trach_nodes(:,2)>max(trach_nodes(:,2))-0.1,:));
%         endpoints = [mean_bot; mean_top];
%         slope_new = (mean_bot(2) - mean_top(2)) / (mean_bot(1) - mean_top(1));
%         carina = [mean_bot, carina(3)];
%         T5     = [mean_body(1), mean_body(2), T5_height];
%         T8     = [mean_body(1), mean_body(2), T8_height];
% 
%         if flags.plot_trachea == 1
%             figure; 
%                 hold on
%                 scatter(trachea_nodes(:,1), trachea_nodes(:,3))
%                 scatter(mean_bot(1), mean_bot(2), "r", "filled")
%                 scatter(mean_top(1), mean_top(2), "r", "filled")
%                 plot(endpoints(:,1), endpoints(:,2),"r")
%                 title(sprintf("Slope of %.2f\n", slope_new))
%         end
% 
%         if sign(slope) ~= sign(slope_new)
%             sign_change = sign_change + 1;
%         end
% 
%         if sign(slope_new) == 1
%             theta = -pi/50 / sign_change;
%         else
%             theta = pi/50 / sign_change;
%         end
% 
%         rotationMatrix = [cos(theta), 0, -sin(theta);...
%                           0,          1,  0;...
%                           sin(theta), 0,  cos(theta)];
% 
%         % Rotate and shift the body
%         nodes  = nodes*rotationMatrix;
%         carina = carina*rotationMatrix;
%         T5     = T5*rotationMatrix;
%         T8     = T8*rotationMatrix;
%     end
% 
%     % Record correct heights from rotation
%     carina_height = carina(3);
%     T5_height     = T5(3);
%     T8_height     = T8(3);

% ----------------------------------------------------------------------- %
%%                          Translate to Origin                           %
% ----------------------------------------------------------------------- %

    if flags.verbose == 1
        fprintf("   Translating to the Origin\n")
    end
    [body_nodes, ~] = Get_Tet_Nodes(nodes, organ_connects{soft_tissue});

    nodes           = nodes         - min(body_nodes);
    sbj_info.carina = carina_height;
    sbj_info.T5     = T5_height;
    sbj_info.T8     = T8_height;

end

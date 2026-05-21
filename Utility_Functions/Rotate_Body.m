function rotated_nodes = Rotate_Body(nodes, body_nodes)
    %{
    Look at the nodes in a torso and rotate them if upside down
    7/8/24 - Kyler Howard

    param: nodes - All nodes in the global mesh
    param: body_nodes - Nodes only from the soft tissue

    return: rotated_nodes - All nodes in the global mesh rotated upside down
    %}

    for theta_d = [0, 180]
        theta_r = theta_d * pi/180;
        % Rotation matrix around z
        % rotationMatrix = [cos(theta_r), -sin(theta_r), 0;...
        %                   sin(theta_r),  cos(theta_r), 0;...
        %                   0,           0,          1];
        % Rotation matrix around y
        rotationMatrix = [cos(theta_r), 0, -sin(theta_r);...
                          0,            1,  0;...
                          sin(theta_r), 0,  cos(theta_r)];
        rotated_body = body_nodes*rotationMatrix;
        % figure; histogram(rotated_body(:,1)); xlabel("x"); title(sprintf("%.2f°", theta_d))
        % figure; histogram(rotated_body(:,3)); xlabel("z"); title(sprintf("%.2f°", theta_d))
        % figure; scatter(rotated_body(:,1), rotated_body(:,3)); xlabel("x"); ylabel("z"); title(sprintf("%.2f°", theta_d))
        fprintf("   %.2f°\n", theta_d)
        fprintf("       Middle: %.2f\n", (min(rotated_body(:,3)) + max(rotated_body(:,3)))/2)
        fprintf("       Mean:   %.2f\n", mean(rotated_body(:,3)))
        z_middle = (min(rotated_body(:,3)) + max(rotated_body(:,3)))/2;
        z_mean   = mean(rotated_body(:,3));

        if z_mean > z_middle
            rotated_nodes      = nodes*rotationMatrix;
            rotated_nodes(:,1) = rotated_nodes(:,1) * -1;
            rotated_nodes(:,3) = rotated_nodes(:,3) + abs(min(rotated_nodes(:,3)));
            fprintf("   Rotated body\n")
            break
        end
    end
end
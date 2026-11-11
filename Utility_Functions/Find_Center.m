function [center, plane] = Find_Center(nodes, height, E, FoB)
    %{
    Find the center of the body and the plane it is on
    9/27/24 - Kyler Howard

    param: nodes  - All nodes on the boundary of the surface
    param: height - The wanted height for the plane
    param: E      - Structure containing the size and type of electrode
    param: FoB    - Flag for if we are looking at the front or back of the body

    return: center - Coordinate of the center
    return: plane  - Plane containing the center
    %}

    % Don't place electrodes higher than the body
    if height > max(nodes(:,3))
        height = max(nodes(:,3));
    end

    % Use whichever one isn't a NaN
    if isnan(E.E_height)
        plane_gap = E.E_rad;
    elseif isnan(E.E_rad)
        plane_gap = E.E_height / 2;
    end

    z_check = (nodes(:,3) > height - plane_gap) & (nodes(:,3) < height + plane_gap);
    plane   = nodes(z_check,:);

    % Find the nodes that are the center of the front/back of the body
    x_mid = (min(plane(:,1)) + max(plane(:,1))) / 2;
    y_min = min(plane(:,2));
    y_max = max(plane(:,2));
    if nargin < 4
        center = [];
    elseif lower(FoB) == "front"
        center = plane(dsearchn(plane, [x_mid, y_max, height]),:);
    elseif lower(FoB) == "back"
        center  = plane(dsearchn(plane, [x_mid, y_min, height]),:);
    end
end
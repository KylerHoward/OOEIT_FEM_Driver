function [centers, normals] = Find_Plot_Normals(connect, nodes)
    normals = zeros(size(connect));
    centers = zeros(size(connect));
    for i = 1:size(connect,1)
        % Extract the corners of the face
        points = nodes(connect(i,:),:);

        % Calculate the incenters and then the center of the face
        side1 = norm(points(2,:) - points(3,:));
        side2 = norm(points(1,:) - points(3,:));
        side3 = norm(points(1,:) - points(2,:));

        centers(i,:) = (side1*points(1,:) + side2*points(2,:) + side3*points(3,:)) / (side1 + side2 + side3);

        % Calculate the unit normal vectors
        normals(i,:) = cross(points(2,:) - points(1,:), points(3,:)-points(1,:));
        normals(i,:) = normals(i,:) / norm(normals(i,:));
    end
end
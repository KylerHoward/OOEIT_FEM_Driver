function radius = Parameratize_Bdry(plane, M, theta)
    %{
    Parameratize a boundary to best place the electrodes, and find the
    radius at a specific point based on a given angle
    9/27/24 - Kyler Howard

    param: plane - A subset of nodes to look through to find centers of electrodes quickly
    param: M     - How many terms to include in the parameratization
    param: theta - What angle (rad) to find the parameratized radius

    return: radius - Parameratized radius at a specific angle
    %}

    % Shift plane to be centered
    center  = (max(plane,[],1) + min(plane,[],1)) / 2;
    c_plane = plane - center;

    % Initialize a radius and theta vector for the plane
    r_vec  = zeros(size(c_plane,1),1);
    th_vec = zeros(size(c_plane,1),1);

    % Loop through every point
    for i = 1:length(c_plane)
        % Calculate the radius from the center to that point
        r_vec(i)  = sqrt(sum((c_plane(i,1:2)).^2));

        % Calculate the angle to that point based on the quadrant
        if c_plane(i,1) > 0 && c_plane(i,2) > 0
            th_vec(i) = atan(c_plane(i,2) / c_plane(i,1));
        elseif c_plane(i,1) < 0 && c_plane(i,2) > 0
            th_vec(i) = pi/2 + atan(abs(c_plane(i,1)) / c_plane(i,2));
        elseif c_plane(i,1) < 0 && c_plane(i,2) < 0
            th_vec(i) = pi + atan(abs(c_plane(i,2)) / abs(c_plane(i,1)));
        elseif c_plane(i,1) > 0 && c_plane(i,2) < 0
            th_vec(i) = 3*pi/2 + atan(c_plane(i,1) / abs(c_plane(i,2)));
        end
    end

    % Create the Q matrix with the sin/cos functions
    Q = NaN(2*M + 1, length(c_plane));
    for i = 0:M
        Q(i+1,:)   = cos(i*th_vec)';
    end
    j = 1;
    for i = (M+2):(2*M+1)
        Q(i,:) = sin(j*th_vec)';
        j = j+1;
    end

    % Solve for the parameters a
    a = (Q*Q')\Q*r_vec;

    % Split the parameters into the initial radius, the cosine
    % coefficients, and the sine coefficients
    r_0 = a(1);
    a_c = a(2:M+1);
    a_s = a(M+2:end);

    % Find the exact value that we are looking for
    radius = r_0;
    for m = 1:M
        radius = radius + a_c(m)*cos(m*theta) + a_s(m)*sin(m*theta);
    end
end
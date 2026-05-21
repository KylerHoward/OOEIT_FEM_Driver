function carina_height = Find_Carina(t_nodes, flags)
    %{
    Determine the height of the carina in the mesh made up of nodes
    8/20/24 - Kyler Howard
    
    param: t_nodes - Tetrahedron node list, n by 3, where n is the number of nodes
    param: flags   - Various flags controlling plotting and other parameters
    
    return: carina_height - Z coordinate of the carina
    %}
    
    % Extract the heights of all nodes
    z_list = unique(t_nodes(:,3));
    
    % Sort the nodes by their z height, and find a cluster radius
    sorted_nodes = sort(t_nodes, 1, "ascend");
    xrange       = range(sorted_nodes(1:round(length(sorted_nodes)/10),1));
    yrange       = range(sorted_nodes(1:round(length(sorted_nodes)/10),2));
    cluster_rad  = max([xrange, yrange]) / 2;
    
    % Find the centroids along the main trachea
    centroids = zeros(0,3);
    i = 1;
    for height = round(max(z_list)):-1:round(max(z_list))-20
        plane = t_nodes(t_nodes(:,3) < height & t_nodes(:,3) > height -1, :);
        centroids(i,:) = mean(plane);
        i = i+1;
    end

    % Use SVD to find a line of best fit
    X_mean = mean(centroids,1);
    X_centered = bsxfun(@minus, centroids, X_mean);
    [~, ~, V] = svd(X_centered, 0);
    direction = V(:,1);
    
    % Loop through at a specific gap
    gap     = 2;
    counter = 1;
    for z = round(max(z_list)):-gap:round(min(z_list))
        plane = t_nodes(t_nodes(:,3) < z & t_nodes(:,3) > z-gap, :);

        % Calculate the parameter t for the given z
        t = (z - X_mean(3)) / direction(3);
        
        % Calculate the corresponding x and y
        x = X_mean(1) + t * direction(1);
        y = X_mean(2) + t * direction(2);

        global_centroid = [x,y,z];
    
        % Calculate the number of clusters
        if size(plane,1) == 1
            clusters = 1;
        else
            clusters = dbscan(plane, cluster_rad, 1, "Distance","squaredeuclidean");
        end
    
        % Find the centroid of each plane and adjust the symbol, color, and
        % name for the gscatter() plotting
        cluster_centroids = zeros(length(unique(clusters)), 3);
        syms   = "";
        colors = 'brgcmy';
        names  = strings(length(unique(clusters)),1);
        for i = unique(clusters).'
            cluster_centroids(i,:) = mean(plane(clusters == i, :), 1);
            syms    = syms + "o";
            names(i) = sprintf("Cluster %d Centroid", i);
        end
        
        % Plot every 10 mm
        if mod(counter*gap, 10) == 0 && flags.plot_slices == 1
            figure()
                hold on
                gscatter(plane(:,1), plane(:,2), clusters, colors(1:length(unique(clusters))), syms)
                gscatter(cluster_centroids(:,1), cluster_centroids(:,2), unique(clusters), colors(1:length(unique(clusters))), "filled")
                scatter(global_centroid(1),    global_centroid(2),   60, "*", "k")
                title(z)
                legend([unique(clusters).', names.', "Global Centroid"], "Location","eastoutside")
        end
        
        % Check if there are 2 clusters and if we've gone past the first 25% of
        % the trachea
    
        % ALTERNATE METHOD: check if the centroid of the previous plane is
        % in-between the two centroids of this plane
        if size(cluster_centroids, 1) == 2% && z > 0.25*max(z_list)
            % Extract the x/y bounds for the two cluster centroids
            x_min = min(cluster_centroids(:,1));
            x_max = max(cluster_centroids(:,1));
            y_min = min(cluster_centroids(:,2));
            y_max = max(cluster_centroids(:,2));
    
            % Check if the global centroid is between either of the clusters
            x_check = (global_centroid(1) > x_min) && (global_centroid(1) < x_max);
            y_check = (global_centroid(2) > y_min) && (global_centroid(2) < y_max);
            if (x_check == 1) || (y_check == 1)
                % Save the carina height and exit the loop
                carina_height = z;
                break
            end
        end
        counter = counter + 1;
    end
    
    fprintf("   The height of the carina is %.2f mm\n", carina_height)
    if flags.plot_trachea == 1
        figure(carina_height)
            hold on
            gscatter(plane(:,1), plane(:,2), clusters, colors(1:length(unique(clusters))), syms)
            gscatter(cluster_centroids(:,1), cluster_centroids(:,2), unique(clusters), colors(1:length(unique(clusters))), "filled")
            scatter(global_centroid(1),    global_centroid(2),   60, "*", "k")
            title(carina_height)
            legend([unique(clusters).', names.', "Global Centroid"], "Location","eastoutside")
    end

end
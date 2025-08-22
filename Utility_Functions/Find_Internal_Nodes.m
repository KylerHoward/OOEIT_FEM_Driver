function inside_indices = Find_Internal_Nodes(t_nodes, flags)
    %{
    Find all internal nodes of a mesh
    8/22/24 - Kyler Howard
    
    param: t_nodes - All nodes in a given tissue
    param: flags   - Various flags controlling plotting and other parameters
    
    return: inside_indices - indicies of internal nodes in a mesh
    %}

    warning('off')

    % How large of a window to look at, at a time (default = 15)
    gap = 10; % mm
    
    % Initialize coords that are inside the boundary
    inside_nodes      = [0, 0, 0];
    inside_node_index = 1;
    for height = floor(min(t_nodes(:,3))) : gap : ceil(max(t_nodes(:,3))-gap)
        % Split into the points we care about
        plane = t_nodes(t_nodes(:,3) > height & t_nodes(:,3) < height + gap, :);
        x     = plane(:,1);
        y     = plane(:,2);
        z     = plane(:,3);

        if isempty(plane)
            continue
        end

        % Find a cluster radius
        cluster_rad  = mean([range(x), range(y)]) / 2;

        % Finding clusters
        clusters = dbscan([x,y], cluster_rad, 1, 'Distance','squaredeuclidean');
        % clusters = dbscan([x,y,z], cluster_rad/2, 5, 'Distance','squaredeuclidean');

        % Find out if a cluster is a min/max
        for group = unique(clusters)'
            if min(y) == min(plane(clusters == group, 2))
                miny_group = group;
            end
            if max(y) == max(plane(clusters == group, 2))
                maxy_group = group;
            end
            if min(x) == min(plane(clusters == group, 1))
                minx_group = group;
            end
            if max(x) == max(plane(clusters == group, 1))
                maxx_group = group;
            end
        end

        % Keep the extreme clusters
        keep_indices   = clusters == miny_group | clusters == maxy_group | clusters == minx_group | clusters == maxx_group;
        outside_points = plane(keep_indices, :);
        inside_points  = plane(~keep_indices,:);

        % Add points to inside list
        for i = 1:size(inside_points,1)
            inside_nodes(inside_node_index,:) = inside_points(i,:);
            inside_node_index = inside_node_index + 1;
        end
            
        % Plot individual slices as looping through
        if flags.plot_slices == 1 && flags.plot_internal
            % Plot the points and the boundary
            figure
                subplot(1,2,1)
                    gscatter(x, y, clusters)
                    xlabel('x')
                    ylabel('y')
                    title(sprintf("Gropus of Points at %d mm", height))
                subplot(1,2,2)
                    hold on
                    scatter(outside_points(:,1),   outside_points(:,2),   12, 'k')
                    scatter(inside_points(:,1), inside_points(:,2), 12, 'r', 'filled')
                    hold off
                    xlabel('x')
                    ylabel('y')
                    title(sprintf("Deleting Red Points at %d mm", height))
        end
    end

    % Plot the original nodes and the internal nodes to be removed
    if flags.plot_internal == 1
        figure
            hold on
            scatter3(t_nodes(:,1), t_nodes(:,2), t_nodes(:,3))
            scatter3(inside_nodes(:,1), inside_nodes(:,2), inside_nodes(:,3),'r','filled')
            xlabel('x')
            ylabel('y')
            zlabel('z')
            title('Found Internal Points')
    end
    
    % Find the rows to delete from nodes, and return them for deletion
    [~, inside_indices, ~] = intersect(t_nodes, inside_nodes, "rows");

end
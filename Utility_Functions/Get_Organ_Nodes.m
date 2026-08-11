function [organ_nodes, organ_faces] = Get_Organ_Nodes(nodes, connectivity, labels, organ_connects, flags)

    % Find the surface shells of each organ
    organ_faces.body  = Get_Unique_Faces(organ_connects.soft_tissue);
    organ_faces.lung  = Get_Unique_Faces(organ_connects.lung);
    organ_faces.heart = Get_Unique_Faces(organ_connects.heart);

    % Account for subjects without a trachea
    if ~isempty(organ_connects.trachea) 
        organ_faces.trachea = Get_Unique_Faces(organ_connects.trachea);
        trachea_intersect   = intersect(sort(organ_faces.body,2),  sort(organ_faces.trachea,2), "rows");
    else
        organ_faces.trachea = [];
        trachea_intersect   = [];
    end

    % Account for subjects without an esophagus
    if ~isempty(organ_connects.esophagus) 
        organ_faces.esophagus = Get_Unique_Faces(organ_connects.esophagus);
        esophagus_intersect   = intersect(sort(organ_faces.body,2),  sort(organ_faces.esophagus,2), "rows");
    else
        organ_faces.esophagus = [];
        esophagus_intersect   = [];
    end

    % Account for subjects without bones
    if ~isempty(organ_connects.bone) 
        organ_faces.bone = Get_Unique_Faces(organ_connects.bone);
        bone_intersect   = intersect(sort(organ_faces.body,2),  sort(organ_faces.bone,2), "rows");
    else
        organ_faces.bone = [];
        bone_intersect   = [];
    end
    
    % Find the intersection of each internal organ and the soft tissue
    lung_intersect  = intersect(sort(organ_faces.body,2),  sort(organ_faces.lung,2),      "rows");
    heart_intersect = intersect(sort(organ_faces.body,2),  sort(organ_faces.heart,2),     "rows");
    
    % Find the difference between the organs to only keep body surface faces
    inner_intersects = vertcat(lung_intersect, trachea_intersect, heart_intersect, esophagus_intersect, bone_intersect);
    surface_faces    = setdiff(sort(organ_faces.body,2), sort(inner_intersects,2), "rows");
    
    % Extract desired organ nodes
    if ~isempty(organ_connects.trachea)
        [organ_nodes.trachea, ~]  = Get_Tet_Nodes(nodes, organ_connects.trachea);
    else
        organ_nodes.trachea = [];
    end
    [organ_nodes.lung, ~]     = Get_Tet_Nodes(nodes, organ_connects.lung);
    [organ_nodes.boundary, ~] = Get_Surface_Nodes(nodes, surface_faces);
    
    if flags.plot_trachea == 1 && ~isempty(organ_connects.trachea)
        figure()
        subplot(1,2,1)
            scatter3(organ_nodes.trachea(:,1), organ_nodes.trachea(:,2), organ_nodes.trachea(:,3), "MarkerEdgeAlpha", 0.2)
            hold on
            scatter3(mean(organ_nodes.trachea(:,1)), mean(organ_nodes.trachea(:,2)), sbj_info.carina, "r", "filled")
            scatter3(organ_nodes.trachea(startsWith(string(organ_nodes.trachea(:,3)), sprintf("%.1f", sbj_info.carina)), 1), organ_nodes.trachea(startsWith(string(organ_nodes.trachea(:,3)), sprintf("%.1f", sbj_info.carina)), 2), sbj_info.carina, "r", "filled")
            % constantplane("z", sbj_info.carina) R2024b
            title(sprintf("Carina: %.2f mm", sbj_info.carina))
            axis equal
        subplot(1,2,2)
            pdeplot3D(nodes.', organ_connects.lung.')
    end
    
    % Find and save the heart node information
    if flags.save_heart_mesh == 1
        [heart_nodes, ~]          = Get_Tet_Nodes(nodes, organ_connects.heart);
        [heart_surface_nodes, ~]  = Get_Surface_Nodes(nodes, organ_faces.heart);
        heart_faces               = organ_faces.heart;
        organ_nodes.heart         = heart_nodes;
        organ_nodes.heart_surface = heart_surface_nodes;
    
        heart_name = sprintf("%s_Heart_Mesh.mat", sbj_name);
        save(fullfile("Heart_Meshes", heart_name), "nodes", "connectivity", "labels", "heart_faces", "heart_nodes", "heart_surface_nodes")
    else
        organ_nodes.heart         = [];
        organ_nodes.heart_surface = [];
    end
end
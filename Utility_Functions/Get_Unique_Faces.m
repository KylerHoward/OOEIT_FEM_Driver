function face_mesh = Get_Unique_Faces(tet_mesh)
    %{
    Determine boundary faces of tetrahedra stored in a mesh
    7/9/24 - Kyler Howard

    param: tet_mesh - tetrahedron index list, m by 4, where m is the number of tetrahedral elements

    return: face_mesh - list of boundary faces, n by 3, where n is the number of boundary faces
    %}
    
    % get all faces
    all_faces = [tet_mesh(:,1), tet_mesh(:,2), tet_mesh(:,3); ...
               tet_mesh(:,1), tet_mesh(:,3), tet_mesh(:,4); ...
               tet_mesh(:,1), tet_mesh(:,4), tet_mesh(:,2); ...
               tet_mesh(:,2), tet_mesh(:,4), tet_mesh(:,3)];
    
    % sort rows so that faces are reordered in ascending order of indices
    sorted_faces = sort(all_faces,2);
    
    % determine uniqueness of faces
    [u,~,idx] = unique(sorted_faces,"rows");

    % determine counts for each unique face
    counts = accumarray(idx(:), 1);

    % extract faces that only occurred once
    sorted_exteriorF = u(counts == 1,:);
    
    % find in original faces so that ordering of indices is correct
    face_mesh = all_faces(ismember(sorted_faces,sorted_exteriorF,"rows"),:);
end

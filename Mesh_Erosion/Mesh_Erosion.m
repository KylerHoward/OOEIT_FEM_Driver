function Mesh_Erosion(filename, filepath, do_prints)
% clear
% clc
% close all
% filename = "EIT212_Segmentation_NoBones.nrrd";
% filepath = "R:\Segmentations\CT_CF_NIH_R01\CTs\EIT212";
% do_prints = 1;

    data   = nrrdinfo(fullfile(filepath, filename));
    volume = nrrdread(fullfile(filepath, filename));
    
    % Set up origin
    origin = data.RawAttributes.spaceorigin;
    origin = str2num(origin(2:end-1));

    step_size = data.RawAttributes.spacedirections;
    step_size = split(step_size, ' ');
    x = str2num(step_size{1}(2:end-1));
    y = str2num(step_size{2}(2:end-1));
    z = str2num(step_size{3}(2:end-1));

    % Create the mesh grid
    X = origin(1) : x(1) : origin(1) + (size(volume,1)-1)*x(1);
    Y = origin(2) : y(2) : origin(2) + (size(volume,2)-1)*y(2);
    Z = origin(3) : z(3) : origin(3) + (size(volume,3)-1)*z(3);

    [Y, X, Z] = meshgrid(Y,X,Z);
    

    % Convert each label into a single binary channel
    if do_prints
        fprintf("   Converting to binary channels\n")
    end
    bwimages = zeros(size(volume,1), size(volume,2), size(volume,3), length(unique(volume)));
    for label = unique(volume)'
        temp = volume;
        temp(temp ~= label) = 0;
        temp = imbinarize(temp);
        bwimages(:,:,:,label+1) = temp;
    end
    clear temp
    
    % Erode/dilate the soft tissue to get rid of any single pixels
    if do_prints
        fprintf("   Eroding/Dilating single pixels of soft tissue\n")
    end
    bwimages(:,:,:,4) = imerode(bwimages(:,:,:,4), strel("cube", 3));
    bwimages(:,:,:,4) = imdilate(bwimages(:,:,:,4), strel("cube", 3));
    
    % Fill the soft tissue. First the top/bottom, then the middle
    if do_prints
        fprintf("   Filling soft tissue\n")
    end
    [~, ~, channels] = ind2sub(size(bwimages(:,:,:,4)), find(bwimages(:,:,:,4) == 1));
    bwimages(:,:,min(channels),4) = imfill(bwimages(:,:,min(channels),4), "holes");
    bwimages(:,:,max(channels),4) = imfill(bwimages(:,:,max(channels),4), "holes");
    bwimages(:,:,:,4) = imfill(bwimages(:,:,:,4), "holes");
    
    % Erode the others
    for label = 2:size(bwimages,4)
        if label == 4
            continue
        end
        if do_prints
            fprintf("   Eroding label %d\n", label)
        end
        bwimages(:,:,:,label) = imerode(bwimages(:,:,:,label), strel("cube", 3));
    end

    % Fill small holes in trachea
    threshold = 25;
    bwimages(:,:,:,3) = bwareaopen(bwimages(:,:,:,3), threshold);
    
    % % Reapply the labels. Had to do 6 before 5 due to the hardcoding of 4
    % if do_prints
        % fprintf("   Reforming an NRRD\n")
    % end
    % fin_vol = bwimages(:,:,:,4);
    % fin_vol(fin_vol == 1) = 3;
    % for label = 2:size(bwimages,4)
    %     if label == 4
    %         continue
    %     end
    %     fin_vol = fin_vol +  bwimages(:,:,:,label);
    %     % Account for label 5 messing everything up
    %     if label == 5
    %         % Set it to something crazy
    %         fin_vol(fin_vol == 4) = 20;
    %     else
    %         fin_vol(fin_vol == 4) = label-1;
    %     end
    % end
    % % Set it back to correct
    % fin_vol(fin_vol == 20) = 4;
    % 
    % if do_prints
    %     fprintf("Visualizing\n")
    % end
    % figure
    % for i = 1:25:size(volume,3)
    %     subplot(1,2,1)
    %         imagesc(volume(:,:,i))
    %         title(sprintf("Original Slice %d", i))
    %         axis off square
    %     subplot(1,2,2)
    %         imagesc(fin_vol(:,:,i))
    %         title(sprintf("Edited Slice %d", i))
    %         axis off square
    %     pause(1)
    % end

    % % Set visualization settings. Delete for looping version
    % color = [0   0   0;...        % Background
    %          246 7   7;...        % Lungs
    %          0   255 0;...        % Trachea
    %          22  46  226;...      % Soft Tissue
    %          222 141 20;...       % Esophagus
    %          222 20  161] ./ 255; % Heart
    % intensity = [-3024 -700 -600 600 700 3071];
    % queryPoints = linspace(min(intensity),max(intensity),256);
    % Cmap = interp1(intensity,color,queryPoints);
    % 
    % % Delete this for looping version
    % vol = volshow(fin_vol,...
    %     RenderingStyle="VolumeRendering",...
    %     Colormap=Cmap);
    
    %%
    
    name_parts = split(filename, '.');
    save_name  = name_parts{1};

    % Loop through each label and save as OBJ
    for label = 2:size(bwimages,4)
        if sum(bwimages(:,:,:,label), 'all') == 0
            continue
        end
        if do_prints
            fprintf("   Saving OBJ %d\n", label-1)
        end

        % Generate a 3D surface from the binary mask
        [faces, vertices] = isosurface(X, Y, Z, bwimages(:,:,:,label), 0.5); % Adjust threshold if needed
        TR = triangulation(faces, vertices);

        % % Define STL file name
        % stlFileName = fullfile(filepath, sprintf('%s_Eroded_Label%d.stl', save_name, label-1));
        % 
        % % Write the STL file
        % stlwrite(TR, stlFileName);

        objFileName = fullfile(filepath, "OBJs", sprintf('%s_Eroded_Label%d.obj', save_name, label-1));
        fid = fopen(objFileName, 'w');

        % Write vertices
        fprintf(fid, '# OBJ file generated from STL\n');
        vertices = TR.Points;
        for i = 1:size(vertices, 1)
            fprintf(fid, 'v %.6f %.6f %.6f\n', vertices(i, 1), vertices(i, 2), vertices(i, 3));
        end

        % Write faces
        faces = TR.ConnectivityList;
        for i = 1:size(faces, 1)
            fprintf(fid, 'f %d %d %d\n', faces(i, 1), faces(i, 2), faces(i, 3));
        end

        % Close the file
        fclose(fid);
    end
end
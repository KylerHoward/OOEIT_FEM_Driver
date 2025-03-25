clc
clear
close all
tic

segment_bones = 0;

% Setup location of Cleaver for Windows PC
if ispc
    cleav_path = "C:\Program Files\Cleaver2\";
    cd(cleav_path)
% Setup location of Cleaver for Apple PC
elseif ismac
    error("I don't know apple file system...")
end

% Open the main folder containing the subjects
atlas_path = uigetdir("Open location of main atlas");
atlas_cont = dir(atlas_path);
run_times  = nan(length(atlas_cont));

% parfor fold_i = 3:length(atlas_cont)
for fold_i = 3:length(atlas_cont)
    % Skip non-folder contents
    if atlas_cont(fold_i).isdir ~= 1
        continue
    end

    sbj_name = atlas_cont(fold_i).name;
    sbj_path = fullfile(atlas_path, sbj_name);
    sbj_cont = dir(sbj_path);

    % Find the segmentation file name
    for sbj_i = 3:length(sbj_cont)
        % Skip folder contents
        if sbj_cont(sbj_i).isdir == 1
            continue
        end

        seg_name = sbj_cont(sbj_i).name;
        if (contains(seg_name, "Segmentation_NoBones.nrrd") == 1) && (segment_bones == 0)
            break
        elseif (contains(seg_name, "Segmentation.nrrd") == 1) && (segment_bones == 1)
            break
        end
    end

    % Create the mesh_name
    msh_name = char(replace(string(seg_name), "Segmentation", "Mesh"));
    % msh_name = replace(string(msh_name), "nrrd", "mat");
    msh_name = msh_name(1:end-5);

    % Create the two full file paths
    seg_path = fullfile(sbj_path, seg_name);
    msh_path = fullfile(sbj_path, "/");

    % % Catch any spaces in file paths
    % if contains(atlas_path, " ")
    %     error("Atlas path/name cannot have spaces")
    % elseif contains(seg_path, " ")
    %     error("Segmentation path/name cannot have spaces")
    % elseif contains(msh_path, " ")
    %     error("mesh path/name cannot have spaces")
    % end

    % Run the Cleaver command line for Windows PC
    if ispc
        tic
        [status, cmdout] = system(sprintf("cleaver-cli -i %s -n %s -f matlab -o %s -e -S", seg_path, msh_name, msh_path));
        if status ~= 0
            error(cmdout)
        end
        run_times(fold_i) = toc;
        fprintf("Created mesh for %s in %.2f minutes\n", sbj_name, run_times(fold_i)/60)

        % Delete the extra info file that was created
        info_name = strcat(msh_name, ".info");
        delete(fullfile(sbj_path, info_name))
    % Run the Cleaver command line for Apple PC
    elseif ismac
        error("I don't know apple file system...")
    end

end
actual_time = toc;

fprintf("\n")
fprintf("Finished all meshes in an total   of %.2f hours\n",   sum(run_times(3:end), 'omitmissing')/3600)
fprintf("Finished all meshes in an average of %.2f minutes\n", mean(run_times(3:end), 'omitmissing')/60)
fprintf("\n")
fprintf("Actual time to compute all meshes is %.2f hours\n", actual_time/3600)
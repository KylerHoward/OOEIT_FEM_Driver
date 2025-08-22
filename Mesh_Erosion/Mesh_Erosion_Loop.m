clear
clc
close all
p = gcp('nocreate'); % Get the current parallel pool without creating a new one
if isempty(p)
    parpool('Processes', 7)
end

do_prints = 0;

% Get contents of server
seg_path     = fullfile("..","..","Segmentations");
seg_contents = dir(seg_path);

sbjcount = 0;
% Loop through all datasets in server
start_time = tic;
for i = 1:length(seg_contents)
    dataset_name = seg_contents(i).name;
    if contains(dataset_name, '.') || seg_contents(i).isdir == 0
        continue
    end

    % Get contents of each dataset
    dataset_path     = fullfile(seg_path, dataset_name, "CTs");
    dataset_contents = dir(dataset_path);

    % Loop through all subjects in each dataset
    parfor ii = 1:length(dataset_contents) 
        clccount = 0;
        sbj_name = dataset_contents(ii).name;
        if contains(sbj_name, '.') || dataset_contents(ii).isdir == 0
            continue
        end

        % Only inspect subjects who have the NRRD file
        nrrd_name = sprintf("%s_Segmentation_NoBones.nrrd", sbj_name);

        nrrd_path = fullfile(dataset_path, sbj_name);
        if isfile(fullfile(nrrd_path, nrrd_name))
            % Make an OBJ folder if it doesn't exist
            if ~isfolder(fullfile(nrrd_path, "OBJs"))
                mkdir(fullfile(nrrd_path, "OBJs"));
            end

            % Only run on subjects who don't have an eroded OBJ already
            obj_name = sprintf("%s_Segmentation_NoBones_Eroded_Label1.obj", sbj_name);
            if ~isfile(fullfile(nrrd_path, "OBJs", obj_name))
                fprintf("%s\n", nrrd_name)
                Mesh_Erosion(nrrd_name, nrrd_path, do_prints)
                clccount = clccount + 1;
                sbjcount = sbjcount + 1;
            end

            % Clear the command line if it gets too long
            if clccount >= 10
                clccount = 0;
                clc
            end
        end
    end
end

stop_time = toc(start_time);
fprintf("It took %.2f minutes to run\n", stop_time/60)
fprintf("An average of %.2f minutes per subject\n", (stop_time/60)/sbjcount)
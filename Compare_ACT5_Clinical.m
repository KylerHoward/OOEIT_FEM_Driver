%{
Compare results of 3D FEM simulation to clinical data
10/8/24 - Kyler Howard

load: Umeas - Measured voltages on the electrodes from FEM
load: current_pat - Current patterns for belt/patch electodes
load: sbj011 - Structure containing voltage and current pattern for subject 011
%}

all_fig = findall(0, "type", "figure");
close(all_fig)
clearvars -except temp msh_path old_msh_path save_path old_save_path first_loc
clc
close all

num_sims       = 2;
do_pauses      = 1;
mean_shift     = 1;
show_clinical  = 1;
show_means     = 0;
% clinical_frame = 71;%[37,44]; % 71 (Subj011)
% clinical_frame = [194,224]; % 207 (Kyler Patch)
clinical_frame = [23619,23664]; % 207 (Kyler Finland)
clinical_frame = 23619;
remove_elec    = 0;
bad_elecs      = [4,5,12];
insp_exp_plot  = 0; % MUST HAVE CLINICAL_FRAMES IN ORDER OF [INSP, EXP] AND ONLY TWO!

% Load files
sbj_path = "C:\Users\kyler\OneDrive\School\Colorado State\Research\Dr. Mueller\FEM\OOEIT_FEM_Driver\Clinical_Data";
% sbj_file = "Sbj42_4x8_patch_circle_2025_02_14_15_10_34_0001"; % Kyler Patch
sbj_file = "Sbj42_25_11_12_15_31_08";                           % Kyler Finland
% sbj_file = "Subj011_2019_09_06_15_29_05_0002";
sbj_data = load(fullfile(sbj_path, sbj_file));
parts    = split(sbj_file,"_");
sbj_name = parts(1);
L = 32;

load_name      = cell(num_sims,1);
load_loc       = cell(num_sims,1);
parts          = cell(num_sims,1);
sim_name       = cell(num_sims,1);
sim_volt       = zeros(L, L-1, num_sims);
sim_CP         = zeros(L, L-1, num_sims);
sim_volt_range = zeros(2, L-1, num_sims);
sim_CP_range   = zeros(2, L-1, num_sims);
for i = 1:num_sims
    % Open the simulated voltage file
    if i == 1
        if exist("first_loc", "var") && ischar(first_loc)
            [load_name{i}, load_loc{i}] = uigetfile(first_loc, sprintf("Open Volt File %d of %d", i, num_sims));
        else
            [load_name{i}, load_loc{i}] = uigetfile("Results/", sprintf("Open Volt File %d of %d", i, num_sims));
            first_loc = load_loc{i};
        end
    else
        [load_name{i}, load_loc{i}] = uigetfile(load_loc{i-1}, sprintf("Open Volt File %d of %d", i, num_sims));
    end

    if load_name{i} == 0
        error("No file selected")
    end

    parts{i} = split(load_name{i}, "_");
    if size(parts{i},1) == 1
        parts{i} = split(load_name{i}, "-");
    end
    
    % Data validation
    while contains(parts{i}{1}, "Volt") ~= 1 && contains(parts{i}{2}, "Volt") ~= 1
        load_name{i} = uigetfile(load_loc{i}, "Wrong File. Open Volt File");
        parts{i} = split(load_name{i}, "_");
        if size(parts{i},1) == 1
            parts{i} = split(load_name{i}, "-");
        end
    end

    % Extract simulation name
    sim_name{i} = string(parts{i}{1}) + ", " + string(parts{i}(contains(parts{i}, 'z')));

    % Load files
    try
        sim_volt(:,:,i) = real(load(fullfile(load_loc{i}, load_name{i})).Umeas);
    catch
        sim_volt(:,:,i) = real(load(fullfile(load_loc{i}, load_name{i})).Umeas_NoNoise);
    end
    try
        sim_CP(:,:,i)   = real(load(fullfile(load_loc{i}, load_name{i})).current_pat);
    catch
        sim_CP(:,:,i)   = real(load(fullfile(load_loc{i}, load_name{i})).cur_pat);
    end

    sim_CP(:,:,i) = sim_CP(:,:,i)*1e6;

    % Set plotting ranges
    sim_volt_range(:,:,i)  = [1.1*min(sim_volt(:,:,i));     1.1*max(sim_volt(:,:,i))];
    sim_CP_range(:,:,i)    = [1.1*min(sim_CP(:,:,i));       1.1*max(sim_CP(:,:,i))];
end

% Remove bad electrodes
if remove_elec == 1
    sim_volt(bad_elecs,:,:) = [];
    sim_CP(bad_elecs,:,:)   = [];
end

%%
% % Split the clinical data
fields       = sort(fieldnames(sbj_data));
idx          = find(strcmp(fields,'frame_voltage'));
DataVol0     = real(sbj_data.(fields{idx}));
DataVol0     = DataVol0 * 1e3; % Convert from V to mV
sbj_Vmulti   = DataVol0(1:L-1, :, :);
[~,~,Slides] = size(sbj_Vmulti);  % Slides is the number of frames
sbj_Vmulti   = permute(sbj_Vmulti,[2,1,3]);  

sbj_volt = zeros(L, L-1, length(clinical_frame));
for i = 1:length(clinical_frame)
    % Extract the frames wanted
    sbj_volt(:,:,i)    = sbj_Vmulti(:,:,clinical_frame(i));
end

sbj_CP = real(sbj_data.(fields{contains(lower(fields), 'cur_pattern')}));
sbj_CP = sbj_CP(:,1:L-1) * 1e6; % Convert from A to uA

% Remove bad electrodes
if remove_elec == 1
    sbj_volt(bad_elecs,:,:) = [];
    sbj_CP(bad_elecs,:,:) = [];
end

% Shift by mean to be 0 centered
if mean_shift == 1
    sbj_volt = sbj_volt - repmat(mean(sbj_volt),size(sbj_volt,1),1);
    sim_volt = sim_volt - repmat(mean(sim_volt),size(sbj_volt,1),1);
end

% Set plotting ranges
sbj_volt_range = [1.1*min(sbj_volt); 1.1*max(sbj_volt)];
sbj_CP_range   = [1.1*min(sbj_CP);   1.1*max(sbj_CP)];

volt_range = zeros(2,L-1,num_sims);
CP_range   = zeros(2,L-1,num_sims);
for i = 1:num_sims
    volt_range(:,:,i) = [min([sim_volt_range(1,:,i); mean(sbj_volt_range(1,:,:), 3)]);...
                         max([sim_volt_range(2,:,i); mean(sbj_volt_range(2,:,:), 3)])];
    CP_range(:,:,i)   = [min([sim_CP_range(1,:,i);   sbj_CP_range(1,1:L-1)]);...
                         max([sim_CP_range(2,:,i);   sbj_CP_range(2,1:L-1)])];
end
volt_range = [min(volt_range(1,:,:), [], 3); max(volt_range(2,:,:), [], 3)];
CP_range   = [min(CP_range(1,:,:), [], 3); max(CP_range(2,:,:), [], 3)];

% Plot each current pattern
figure
for CP = [1:L-1, 1]
% for CP = [1, 18]
    s1 = subplot(1,2,1);
        cla(s1)
        hold on
        if show_clinical == 1
            for i = 1:length(clinical_frame)
                p1 = plot(sbj_volt(:,CP,i) - mean(sbj_volt(:,CP,i)));
                color1 = get(p1, 'Color');
                if show_means == 1
                    plot(mean(sbj_volt(:,CP,i)) * ones(L-1,1),'Color',color1,'LineStyle',':')
                end
            end
        end
        for i = 1:num_sims
            p2 = plot(sim_volt(:,CP,i) - mean(sim_volt(:,CP,i)));
            color2 = get(p2, 'Color');
            if show_means == 1
                plot(mean(sim_volt(:,CP,i)) * ones(L-1,1),'Color',color2,'LineStyle',':')
            end
        end
        title(sprintf("Voltages for Current Pattern %d", CP))
        ylim(volt_range(:,CP))
        ylabel("mV")
        xlabel("Electrode")
        names = [];
        if show_clinical == 1 && show_means == 1
            for i = 1:length(clinical_frame)
                names = [names, sprintf("%s, Frame %d", sbj_name, clinical_frame(i)), sprintf("%s Mean Value, Frame %d", sbj_name, clinical_frame(i))];
            end
        elseif show_clinical == 1 && show_means == 0
            for i = 1:length(clinical_frame)
                names = [names, sprintf("%s, Frame %d", sbj_name, clinical_frame(i))];
            end
        else
        end
        if show_means == 1
            for i = 1:num_sims
                names = [names, sim_name{i}, sprintf("%s Mean Value", sim_name{i})];
            end
        else
            for i = 1:num_sims
                names = [names, sim_name{i}];
            end
        end
        legend(names, 'Location','southoutside')
        axis square
    s2 = subplot(1,2,2);
        cla(s2)
        hold on
        if show_clinical == 1
            plot(sbj_CP(:,CP),'Color',color1)
        end
        for i = 1:num_sims
            plot(sim_CP(:,CP,i),'Color',color2);
        end
        title(sprintf("Current Pattern %d", CP))
        ylim(CP_range(:,CP))
        ylabel("µA")
        xlabel("Electrode")
        if show_clinical == 1
            names = sbj_name;
        else
            names = [];
        end
        for i = 1:num_sims
            names = [names, sim_name{i}];
        end
        legend(names, 'Location','southoutside')
        axis square
    
    if do_pauses == 1
        pause
    end
end

%% Calculate scale difference
if remove_elec == 1
    L = L - length(bad_elecs);
end

volt_scale = zeros(L,num_sims);
CP_scale   = zeros(L,num_sims);
error      = zeros(num_sims,1);
el_error   = zeros(L, num_sims);
for i = 1:num_sims
    for l = 1:L
        volt_scale(l,i) = mean(rmoutliers(mean(sbj_volt(l,:,:),3) ./ sim_volt(l,:,i)));
        CP_scale(l,i)   = mean(rmoutliers(sim_CP(l,:,i)   ./ sbj_CP(l,:)),"all","omitmissing");

        el_error(l,i)   = norm(squeeze(sim_volt(l,:,i) - sbj_volt(l,:,:)), inf) / norm(squeeze(sbj_volt(l,:,:)), inf) * 100;
    end

    vect_sim_volt = sim_volt(:,:,i);
    vect_sbj_volt = sbj_volt(:,:,:);
    error(i) = norm(vect_sim_volt(:) - vect_sbj_volt(:), inf)/norm(vect_sbj_volt(:), inf) * 100;
end

% for i = 1:num_sims
%     figure()
%     hold on
%     bar(1:L, volt_scale(:,i));
%     plot(1:L, mean(volt_scale(:,i))*ones(1,L),'r')
%     xlabel("Electrode")
%     title("Voltage Scale")
%     legend(sim_name{i}, sprintf("%s Mean Scale", sim_name{i}))
% end

for i = 1:num_sims
    fprintf("The average scale difference is %.2f & %.2f for voltage and CP respectively for %s\n", mean(volt_scale(:,i)), mean(CP_scale(:,i)), load_name{i})
    fprintf("The mean electrode  Linf voltage error is %.2f%% for %s\n", mean(el_error(:,i)), load_name{i})
    fprintf("The overal relative Linf voltage error is %.2f%% for %s\n\n", error(i), load_name{i})
end

if insp_exp_plot == 1
    sim_diff = sim_volt(:,:,1) - sim_volt(:,:,2);
    sbj_diff = sbj_volt(:,:,1) - sbj_volt(:,:,2);

    fprintf("Norm difference of simulations is %.2f mV\n", norm(sim_diff))
    fprintf("Norm difference of clinical    is %.2f mV\n", norm(sbj_diff))
end

% find the minimums
[min_el_error_val, min_el_error_idx] = min(mean(el_error));
[min_error_val,    min_error_idx]    = min(error);
fprintf("%s\n", repmat('-',[1,100]))
fprintf("Minimum mean electrode  Linf voltage error is %.2f%% for %s\n", min_el_error_val, load_name{min_el_error_idx})
fprintf("Minimum overal relative Linf voltage error is %.2f%% for %s\n", min_error_val,    load_name{min_error_idx})
fprintf("%s\n", repmat('-',[1,100]))

for i = 1:num_sims
    figure
        hold on
        bar(1:size(el_error,1), el_error(:,i))
        yline(mean(el_error(:,i)), 'r', 'LineWidth',2)
        xlabel("Electrode")
        ylabel("Linf Relative Error")
        title(sprintf("Electrode Error for %s", sim_name{i}))
end

return
% electrode_errors is a 32xn array of the el_error on each subject
mean_vals = mean(el_error,2);
std_vals  = std(el_error, [], 2);

figure; 
    hold on
    plot(1:32, mean_vals, 'r-', 'LineWidth', 2)
    plot(1:32, mean_vals-std_vals, 'b--', 'LineWidth', 1.5)
    plot(1:32, mean_vals+std_vals, 'b--', 'LineWidth', 1.5)
    xlabel("Electrode", "FontSize", 15) 
    ylabel("L-Infinity Relative Error (%)", "FontSize", 15)
    legend(["Mean", "-1 SD", "+1 SD"], "Location", "northeast", "FontSize", 15)

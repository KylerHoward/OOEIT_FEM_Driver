%{
Compare results of 3D FEM simulation to clinical data
10/8/24 - Kyler Howard

load: Umeas - Measured voltages on the electrodes from FEM
load: current_pat - Current patterns for belt/patch electodes
load: subj011 - Structure containing voltage and current pattern for subject 011
%}

clear
clc
close all

num_sims       = 2;
do_pauses      = 0;
mean_shift     = 1;
show_clinical  = 1;
show_means     = 0;
% clinical_frame = 71;%[37,44]; % 71 (Subj011)
clinical_frame = [194,224]; % 207 (Kyler)

% Load files
subj_path = "C:\Users\kyler\OneDrive\School\Colorado State\Research\Dr. Mueller\FEM\OOEIT_FEM_Driver\Clinical_Data";
subj_file = "Sbj42_4x8_patch_circle_2025_02_14_15_10_34_0001";
% subj_file = "Subj011_2019_09_06_15_29_05_0002";
subj_data = load(fullfile(subj_path, subj_file));
parts     = split(subj_file,"_");
subj_name = parts(1);
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
        [load_name{i}, load_loc{i}] = uigetfile("Results/", "Open Volt File");
    else
        [load_name{i}, load_loc{i}] = uigetfile(load_loc{i-1}, "Open Volt File");
    end

    if load_name{i} == 0
        error("No file selected")
    end

    parts{i} = split(load_name{i}, "_");
    
    % Data validation
    while contains(parts{i}{1}, "Volt") ~= 1
        load_name{i} = uigetfile(load_loc{i}, "Wrong File. Open Volt File");
        parts{i} = split(load_name{i}, "_");
    end

    % Extract simulation name
    sim_name{i} = string(parts{i}{2}(1:end-4));

    % Load files
    sim_volt(:,:,i) = real(load(fullfile(load_loc{i}, load_name{i})).Umeas);
    try
        sim_CP(:,:,i)   = real(load(fullfile(load_loc{i}, load_name{i})).current_pat);
    catch
        sim_CP(:,:,i)   = real(load(fullfile(load_loc{i}, load_name{i})).current_pattern);
    end

    % Set plotting ranges
    sim_volt_range(:,:,i)  = [1.1*min(sim_volt(:,:,i));     1.1*max(sim_volt(:,:,i))];
    sim_CP_range(:,:,i)    = [1.1*min(sim_CP(:,:,i));       1.1*max(sim_CP(:,:,i))];
end

%%
% % Split the clinical data
fields       = sort(fieldnames(subj_data));
DataVol0     = real(subj_data.(fields{1}));
[Slides,n]   = size(DataVol0);  % Slides is the number of frames, n the number of measurements
subj_Vmulti  = DataVol0.';      % Change to num_mesh_elts by number of frames
subj_Vmulti  = reshape(subj_Vmulti,L,L-1,Slides);
subj_Vscale  = real(subj_data.(fields{contains(lower(fields), 'vscale')}))';
subj_Vscale  = subj_Vscale(1:L);

% Scale the voltages
subj_Vscale(L+1:end) = [];
subj_Vscale = repmat(subj_Vscale,1,L-1);
for frame = 1:Slides
    subj_Vmulti(:,:,frame) = squeeze(subj_Vmulti(:,:,frame)).*subj_Vscale;
end
    

subj_volt = zeros(L, L-1, length(clinical_frame));
for i = 1:length(clinical_frame)
    % Extract the frames wanted
    subj_volt(:,:,i)    = subj_Vmulti(:,:,clinical_frame(i));
end

subj_CPscale = real(subj_data.(fields{contains(lower(fields), 'iscale')}));
subj_CP      = real(subj_data.(fields{contains(lower(fields), 'pattern')}));
subj_CP      = mean(subj_CPscale)* subj_CP';

% Shift by mean to be 0 centered
if mean_shift == 1
    subj_volt = subj_volt - repmat(mean(subj_volt),32,1);
    sim_volt     = sim_volt     - repmat(mean(sim_volt),32,1);
end

% Set plotting ranges
subj_volt_range = [1.1*min(subj_volt); 1.1*max(subj_volt)];
subj_CP_range   = [1.1*min(subj_CP);   1.1*max(subj_CP)];

volt_range = zeros(2,L-1,num_sims);
CP_range   = zeros(2,L-1,num_sims);
for i = 1:num_sims
    volt_range(:,:,i) = [min([sim_volt_range(1,:,i); mean(subj_volt_range(1,:,:), 3)]);...
                         max([sim_volt_range(2,:,i); mean(subj_volt_range(2,:,:), 3)])];
    CP_range(:,:,i)   = [min([sim_CP_range(1,:,i);   subj_CP_range(1,1:L-1)]);...
                         max([sim_CP_range(2,:,i);   subj_CP_range(2,1:L-1)])];
end
volt_range = [min(volt_range(1,:,:), [], 3); max(volt_range(2,:,:), [], 3)];
CP_range   = [min(CP_range(1,:,:), [], 3); max(CP_range(2,:,:), [], 3)];

% Plot each current pattern for frame 71
figure
for CP = [1:L-1, 1]
% for CP = [1, 18]
    s1 = subplot(1,2,1);
        cla(s1)
        hold on
        if show_clinical == 1
            for i = 1:length(clinical_frame)
                p1 = plot(subj_volt(:,CP,i) - mean(subj_volt(:,CP,i)));
                color1 = get(p1, 'Color');
                if show_means == 1
                    plot(mean(subj_volt(:,CP,i)) * ones(L-1,1),'Color',color1,'LineStyle',':')
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
                names = [names, sprintf("%s, Frame %d", subj_name, clinical_frame(i)), sprintf("%s Mean Value, Frame %d", subj_name, clinical_frame(i))];
            end
        elseif show_clinical == 1 && show_means == 0
            for i = 1:length(clinical_frame)
                names = [names, sprintf("%s, Frame %d", subj_name, clinical_frame(i))];
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
            plot(subj_CP(:,CP),'Color',color1)
        end
        for i = 1:num_sims
            plot(sim_CP(:,CP,i),'Color',color2);
        end
        title(sprintf("Current Pattern %d", CP))
        ylim(CP_range(:,CP))
        ylabel("µA")
        xlabel("Electrode")
        if show_clinical == 1
            names = subj_name;
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
volt_scale = zeros(L,num_sims);
CP_scale   = zeros(L,num_sims);
for i = 1:num_sims
    for l = 1:L
        volt_scale(l,i) = mean(rmoutliers(mean(subj_volt(l,:,:),3) ./ sim_volt(l,:,i)));
        CP_scale(l,i)   = mean(rmoutliers(sim_CP(l,:,i)   ./ subj_CP(l,1:L-1)),"all","omitmissing");
    end
end

for i = 1:num_sims
    figure()
    hold on
    bar(1:L, volt_scale(:,i));
    plot(1:L, mean(volt_scale(:,i))*ones(1,L),'r')
    xlabel("Electrode")
    title("Voltage Scale")
    legend(sim_name{i}, sprintf("%s Mean Scale", sim_name{i}))
end

for i = 1:num_sims
    fprintf("The average scale difference is %.2f & %.2f for voltage and CP respectively for %s\n", mean(volt_scale(:,i)), mean(CP_scale(:,i)), sim_name{i})
end

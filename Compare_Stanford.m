%{
Compare results of 3D FEM simulation to clinical data
10/8/24 - Kyler Howard

load: Umeas - Measured voltages on the electrodes from FEM
load: current_pat - Current patterns for belt/patch electodes
load: subj011 - Structure containing voltage and current pattern for subject 011
%}

clear
% clc
close all

num_sims       = 2;
do_pauses      = 0;
mean_shift     = 1;
show_clinical  = 1;
show_means     = 0;
clinical_frame = [111,117]; % Sbj005
% clinical_frame = [1305,1310]; % Sbj002
% clinical_frame = [1280,1287]; % Sbj003
% clinical_frame = [1217,1221]; % Sbj004
% clinical_frame = 1287;
remove_elec    = 0;
bad_elecs      = [4,5,12];
insp_exp_plot  = 1; % MUST HAVE CLINICAL_FRAMES IN ORDER OF [INSP, EXP] AND ONLY TWO!

% Load files
sbj_path = "C:\Users\kyler\OneDrive\School\Colorado State\Research\Dr. Mueller\FEM\OOEIT_FEM_Driver\Clinical_Data";
% sbj_file = "Sbj42_4x8_patch_circle_2025_02_14_15_10_34_0001";
% sbj_file = "Subj011_2019_09_06_15_29_05_0002";
sbj_file = "ETT_005_circular2x16_11_34_40_0001";
% sbj_file = "ETT_002_circular2x16_12_19_48_0000";
% sbj_file = "ETT_003_circular4x8_11_41_35_0003";
% sbj_file = "ETT_004_circular4x8_13_13_57_0006";
sbj_data = load(fullfile(sbj_path, sbj_file));
parts    = split(sbj_file,"_");
sbj_name = parts(2);
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
    sim_name{i} = string(parts{i}{1});

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
fields      = sort(fieldnames(sbj_data));
DataVol0    = real(sbj_data.(fields{1}));
[Slides,n]  = size(DataVol0);  % Slides is the number of frames, n the number of measurements
sbj_Vmulti  = DataVol0.';      % Change to num_mesh_elts by number of frames
sbj_Vmulti  = reshape(sbj_Vmulti,L,L-1,Slides);
sbj_Vscale  = real(sbj_data.(fields{contains(lower(fields), 'vscale')}))';
sbj_Vscale  = sbj_Vscale(1:L);

% Scale the voltages
sbj_Vscale(L+1:end) = [];
sbj_Vscale = repmat(sbj_Vscale,1,L-1);
for frame = 1:Slides
    sbj_Vmulti(:,:,frame) = squeeze(sbj_Vmulti(:,:,frame)).*sbj_Vscale;
end
    

sbj_volt = zeros(L, L-1, length(clinical_frame));
for i = 1:length(clinical_frame)
    % Extract the frames wanted
    sbj_volt(:,:,i)    = sbj_Vmulti(:,:,clinical_frame(i));
end

sbj_CPscale = real(sbj_data.(fields{contains(lower(fields), 'iscale')}));
sbj_CP      = real(sbj_data.(fields{contains(lower(fields), 'pattern')}));
sbj_CP      = mean(sbj_CPscale)* sbj_CP';

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
figure()
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
        if insp_exp_plot == 1
            names = [names, sprintf("%s Inspiration", sbj_name), sprintf("%s Expiration", sbj_name)];
        else
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
    elseif do_pauses < 1 & do_pauses > 0
        pause(0.75)
    end
end

%% Calculate scale difference
if remove_elec == 1
    L = L - length(bad_elecs);
end
volt_scale = zeros(L,num_sims);
CP_scale   = zeros(L,num_sims);
error      = zeros(num_sims,1);
for i = 1:num_sims
    for l = 1:L
        volt_scale(l,i) = mean(rmoutliers(mean(sbj_volt(l,:,:),3) ./ sim_volt(l,:,i)));
        CP_scale(l,i)   = mean(rmoutliers(sim_CP(l,:,i)   ./ sbj_CP(l,:)),"all","omitmissing");
    end

    error(i) = norm(sim_volt(:,:,i) - sbj_volt(:,:,1))/norm(sbj_volt(:,:,1)) * 100;
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
    fprintf("The relative voltage error is %.2f%% for %s\n\n", error(i), load_name{i})
end

if insp_exp_plot == 1
    sim_diff = sim_volt(:,:,1) - sim_volt(:,:,2);
    sbj_diff = sbj_volt(:,:,1) - sbj_volt(:,:,2);

    fprintf("Norm difference of simulations is %.2f mV\n", norm(sim_diff))
    fprintf("Norm difference of clinical    is %.2f mV\n", norm(sbj_diff))
end


% Analyze the data from patient 104 - piece all of the data sets together
% and plot the PCA
clear all;
close all
clc;


% First seven sets are with ventilator, supine AC, VC+, 40%, 5 PEEP
%% Loading HST Dataset
load Subj011_2019_09_06_15_29_05_0002.mat ;
DataVol0=real(Subj011_2019_09_06_15_29_05_0002); 
[m,n]=size(DataVol0);  % m is the number of frames, n the number of measurements
Slides=m;                  % Slides = total frames of the dataset
Vmulti_t = DataVol0.';  %Change to num_mesh_elts by number of frames
Vmulti = reshape(Vmulti_t,32,31,Slides);

Vmulti_CP1=squeeze(real(Vmulti(:,1,1)));
figure
plot(Vmulti_CP1)
title('Voltages for CP 1')

%% Define applied currents - they will be the same every time
current_pat = Subj011_2019_09_06_15_29_05_0002_Pattern';

% Plot voltages, looping through CPs
% Use frame 71
figure
for ii = 1:31
    subplot(1,2,1)
plot(real(Vmulti(:,ii,71)));
title(['Voltage for CP ', num2str(ii)])
axis square
subplot(1,2,2)
plot(current_pat(:,ii))
title(['CP ', num2str(ii)])
axis square
pause
end

% [Slides,num_mesh_elts]=size(DataVol0);                  % Slides = total frames of the dataset
% [first_comp,ref_frame] = findRefFrameTest(real(DataVol0),Slides,32);

% Compute power waveform
for ii=1:Slides
  Vframe= squeeze(real(Vmulti(:,:,ii)));
  power_wavef_mx(:,:,ii)=current_pat'*Vframe; % should be size 31 by 31 by num frames
end
figure
power_wavef = squeeze(power_wavef_mx(1,1,:));
plot(power_wavef)
title('Power Waveform First data set')
return

%% Loading HST Dataset
load Subj011_2019_09_06_15_34_14_0003.mat ;
DataVol1=real(Subj011_2019_09_06_15_34_14_0003); 
[m,n]=size(DataVol1);  % m is the number of frames, n the number of measurements
Slides=m;                  % Slides = total frames of the dataset
Vmulti_t = DataVol1.';  %Change to num_mesh_elts by number of frames
Vmulti = reshape(Vmulti_t,32,31,Slides);

Vmulti_CP1=squeeze(real(Vmulti(:,1,:)));
figure
mesh(Vmulti_CP1)
title('Voltages for CP 1')

% % Plot voltages, looping through CPs
% % Use frame 200
% figure
% for ii = 1:31
%     subplot(1,2,1)
% plot(real(Vmulti(:,ii,200)));
% title(['Voltage for CP ', num2str(ii)])
% axis square
% subplot(1,2,2)
% plot(current_pat(:,ii))
% title(['CP ', num2str(ii)])
% axis square
% pause
% end

% Compute power waveform
for ii=1:Slides
  Vframe= squeeze(real(Vmulti(:,:,ii)));
  power_wavef_mx(:,:,ii)=current_pat'*Vframe; % should be size 31 by 31 by num frames
end
figure
power_wavef = squeeze(power_wavef_mx(1,1,:));
plot(power_wavef)
title('Power Waveform Second data set')

%return


%% Loading HST Dataset
load Subj011_2019_09_06_15_45_02_0004.mat ;
DataVol3=real(Subj011_2019_09_06_15_45_02_0004); 
[m,n]=size(DataVol3);  % m is the number of frames, n the number of measurements
Slides=m;                  % Slides = total frames of the dataset
Vmulti_t = DataVol3.';  %Change to num_mesh_elts by number of frames
Vmulti = reshape(Vmulti_t,32,31,Slides);

Vmulti_CP1=squeeze(real(Vmulti(:,1,:)));
figure
mesh(Vmulti_CP1)
title('Voltages for CP 1')

% Plot voltages, looping through CPs
% Use frame 200
% figure
% for ii = 1:31
%     subplot(1,2,1)
% plot(real(Vmulti(:,ii,200)));
% title(['Voltage for CP ', num2str(ii)])
% axis square
% subplot(1,2,2)
% plot(current_pat(:,ii))
% title(['CP ', num2str(ii)])
% axis square
% pause
% end

% Compute power waveform
for ii=1:Slides
  Vframe= squeeze(real(Vmulti(:,:,ii)));
  power_wavef_mx(:,:,ii)=current_pat'*Vframe; % should be size 31 by 31 by num frames
end
figure
power_wavef = squeeze(power_wavef_mx(1,1,:));
plot(power_wavef)
title('Power Waveform Fourth data set')

%return

%% Loading HST Dataset
load Subj011_2019_09_06_15_52_28_0005.mat ;
DataVol4=real(Subj011_2019_09_06_15_52_28_0005);
[m,n]=size(DataVol4);  % m is the number of frames, n the number of measurements
Slides=m;                  % Slides = total frames of the dataset
Vmulti_t = DataVol4.';  %Change to num_mesh_elts by number of frames
Vmulti = reshape(Vmulti_t,32,31,Slides);

Vmulti_CP1=squeeze(real(Vmulti(:,1,:)));
figure
mesh(Vmulti_CP1)
title('Voltages for CP 1')

% Plot voltages, looping through CPs
% Use frame 200
figure
for ii = 1:31
    subplot(1,2,1)
plot(real(Vmulti(:,ii,200)));
title(['Voltage for CP ', num2str(ii)])
axis square
subplot(1,2,2)
plot(current_pat(:,ii))
title(['CP ', num2str(ii)])
axis square
pause
end

% Compute power waveform
for ii=1:Slides
  Vframe= squeeze(real(Vmulti(:,:,ii)));
  power_wavef_mx(:,:,ii)=current_pat'*Vframe; % should be size 31 by 31 by num frames
end
figure
power_wavef = squeeze(power_wavef_mx(1,1,:));
plot(power_wavef)
title('Power Waveform Fifth data set')
return


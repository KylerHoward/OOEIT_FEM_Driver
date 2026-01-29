clearvars -except filepath filename
clc
close all

if exist('filepath', 'var') && ischar(filepath)
    [filename, filepath] = uigetfile(filepath, "Open Voltage File to Inspect");
else
    [filename, filepath] = uigetfile("Open Voltage File to Inspect");
end

loadfile = load(fullfile(filepath, filename));
fields   = fieldnames(loadfile);

if contains(fields{1}, "Scaled") == 1 % must be GE data
    % Load and scale the voltage
    fields  = sort(fields);
    voltage = eval(sprintf('loadfile.%s', fields{1}));
    Vscale  = eval(sprintf('loadfile.%s_VScale', fields{1}));
    Vscale  = mean(Vscale(1:32));
    voltage = voltage .* Vscale;
    voltage = voltage.';
    voltage = reshape(voltage, [32,31,size(voltage,2)]);
    voltage = permute(voltage, [2,1,3]);
    voltage = real(voltage);

    % Load and scale the current
    cur_pattern = eval(sprintf('loadfile.%s_Pattern', fields{1}));
    Iscale      = eval(sprintf('loadfile.%s_IScale',  fields{1}));
    cur_pattern = cur_pattern .* mean(Iscale(1:32));
    cur_pattern = cur_pattern'/1000; % Scale it for the later scaling back
else % must be ACT5 data
    voltage = loadfile.frame_voltage;
    voltage = squeeze(real(voltage(1:31,:,:))) * 1000;

    cur_pattern = loadfile.cur_pattern;
end
%%
J = real(cur_pattern(:,1:31));  % These are in amps and they already include the amplitude
cpnormvec = zeros(1,31);
J=J*1000;  % Convert to mA
for kk = 1:31
 cpnormvec(kk) = norm(J(:,kk),2);
 J(:,kk) = J(:,kk)/cpnormvec(kk);
end

power_wavef_mx = zeros(31,31,size(voltage,3));
% Compute power waveform
for ii=1:size(voltage,3)
  Vframe= squeeze(real(voltage(:,:,ii)));
  power_wavef_mx(:,:,ii)=J'*Vframe.'; % should be size 31 by 31 by num frames
end

figure(WindowState='maximized')
    power_wavef = squeeze(power_wavef_mx(1,1,:));
    plot(power_wavef)
    % plot(squeeze(voltage(1,1,:)))
    xlabel("Frame")
    ylabel("mV")
    title(sprintf("Power Waveform of %s", filename), 'Interpreter', 'none')
%%
insp_frame = input("Type the frame number for peak inspiration: ");
exp_frame  = input("Type the frame number for peak expiration:  ");

insp_volt = squeeze(voltage(:,:,insp_frame));
exp_volt  = squeeze(voltage(:,:,exp_frame));

% Shift to 0 DC4
insp_volt = insp_volt - mean(insp_volt,2);
exp_volt  = exp_volt - mean(exp_volt,2);

% figure
%     hold on
%     plot(insp_volt(1,:))
%     plot(exp_volt(1,:))
%     xlabel("Electrode")
%     ylabel("mV")
%     title("Voltage for Current Pattern 1")
%     legend("Inspiration", "Expiration")

volt_diff = insp_volt - exp_volt;
fprintf("Norm difference is %.2f mV\n", norm(volt_diff))
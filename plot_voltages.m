clearvars -except filepath
clc
% close all
    
if exist('filepath', 'var') && ischar(filepath)
    [filename, filepath] = uigetfile(filepath, "Open Voltage File to Inspect");
else
    [filename, filepath] = uigetfile("Open Voltage File to Inspect");
end

fprintf("Inspecting %s\n", filename)
load(fullfile(filepath, filename))  

% Truncate the number of frames
try % Try ACT5 loading
    frame_voltage = squeeze(frame_voltage(1:31,:,:));  % CP, electrode, frame
catch 
    try % GE loading. Must multiply by Vscale/Iscale
        frame_voltage = eval(squeeze(filename(1:end-4)));
        Vscale        = eval(sprintf('%s_VScale', filename(1:end-4)));
        Vscale        = mean(Vscale(1:32));
        frame_voltage = frame_voltage .* Vscale;
        frame_voltage = reshape(frame_voltage, [size(frame_voltage,1),32,31]);
        frame_voltage = permute(frame_voltage, [3,2,1]);
    
        cur_pattern = eval(sprintf('%s_Pattern', filename(1:end-4)));
        Iscale      = eval(sprintf('%s_IScale',  filename(1:end-4)));
        cur_pattern = cur_pattern .* mean(Iscale(1:32));
        cur_pattern = cur_pattern'/1000; % Scale it for the later scaling back
        current_amp = max(real(cur_pattern),[],'all');
    catch % FEM Data
        frame_voltage = permute(Umeas_NoNoise, [2,1,3]);
        cur_pattern   = cur_pat;
        current_amp   = max(real(cur_pattern),[],'all');
    end
end
[K,L,Slides] = size(frame_voltage);
try
    minutes = floor(Slides  / system_frame_rate / 60);
    seconds = floor((Slides / system_frame_rate / 60 - minutes)*60);
catch
    try % ACT 5
        minutes = floor(Slides  / frame_rate / 60);
        seconds = floor((Slides / frame_rate / 60 - minutes)*60);   
    catch 
        try % GE
            system_frame_rate = eval(sprintf('%s_FPS', filename(1:end-4)));
            minutes = floor(Slides  / system_frame_rate / 60);
            seconds = floor((Slides / system_frame_rate / 60 - minutes)*60);
        catch % FEM
            system_frame_rate = flags.fps;
            minutes = floor(Slides  / system_frame_rate / 60);
            seconds = floor((Slides / system_frame_rate / 60 - minutes)*60);
        end
    end
end
    
fprintf("The voltage is %d by %d by %d\n", K, L, Slides)
fprintf("   That is %d minutes and %d seconds\n", minutes, seconds)

% figure
% plot(real(squeeze(frame_voltage(1,1,:))))
% return

%=========================== Define Parameters ============================
V_nonhom = squeeze(real(frame_voltage(1:31,:,10)))*1000;     % mV
%V_hom = squeeze(real(frame_voltage(1:31,:,20)))*1000;    % mV
% Construct voltage matrix. Rows = electrode num , cols= current pattern
V = V_nonhom';           % Reshape voltages into matrix 
%Vref = V_hom';           % Reshape voltages into matrix


CurrAmp = current_amp*1000;             % Current amplitude (mA)
fprintf("The current applitude is %.2f mA\n", CurrAmp)
  J = real(cur_pattern(:,1:31));  % These are in amps and they already include the amplitude
  J=J*1000;  % Convert to mA
  for kk=1:31
     cpnormvec(kk) = norm(J(:,kk),2);
     J(:,kk) = J(:,kk)/cpnormvec(kk);
  end
% end


% Normalize the entries of V so that the voltages sum to zero in each col.
S = sum(V)/L; 
adjust = zeros(L,L-1);
adjust(1:L-1,:) = meshgrid(S);
adjust(L,:) = S;
V = V - adjust;

% Scale columns of V to match L2-normalization of the C.P. matrix J
% if (CP_flag==0)
%   V = V * sqrt(2/L)/CurrAmp;
%   V(:,L/2) = V(:,L/2) * sqrt(1/2); % The L/2 col. gets different treatment
% else
     for kk=1:31
        V(:,kk) = V(:,kk)/cpnormvec(kk); 
     end
% end

% % Normalize the entries of Vref so that the voltages sum to zero in each col.
% S = sum(Vref)/L; 
% adjust = zeros(L,L-1);
% adjust(1:L-1,:) = meshgrid(S);
% adjust(L,:) = S;
% Vref = Vref - adjust;
% 
% % Scale columns of V to match L2-normalization of the C.P. matrix J
% % if (CP_flag==0)
% %   Vref = Vref * sqrt(2/L)/CurrAmp;
% %   Vref(:,L/2) = V(:,L/2) * sqrt(1/2); % The L/2 col. gets different treatment
% % else
%      for kk=1:31
%         Vref(:,kk) = Vref(:,kk)/cpnormvec(kk);
%      end
% % end



R = V' * J;                     % R is the ND map 
Lambda = inv(R);                % Lambda is the DN map, size L-1 x L-1

% Rref = Vref' * J;                     % R is the ND map 
% homLam = inv(Rref);                % Lambda is the DN map, size L-1 x L-1

% figure
% mesh(real(R))
% figure
% % plot them
% for ii = 1:31
%     plot(J(:,ii));
%     title('Normalized CPs')
%     pause 
% end
% figure
% % plot them
% for ii = 1:31
%     plot(V(:,ii));
%     title('Normalized Voltages')  % Check old code and see if these sizes look right 
% end
%return

power_wavef_mx = zeros(31,31,Slides);
% Compute power waveform
for ii=1:Slides
  Vframe= squeeze(real(frame_voltage(:,:,ii)));
  power_wavef_mx(:,:,ii)=J'*Vframe.'; % should be size 31 by 31 by num frames
end
figure
    hold on
    power_wavef = squeeze(power_wavef_mx(1,1,:));
    plot(power_wavef)
    title(sprintf("Power Waveform of %s", filename), 'Interpreter', 'none')

    if exist('annotations', 'var')
        fprintf("Annotations:\n")
        for i = 1:length(annotations)
            if isempty(annotations(i).frame) == 0
                xline(annotations(i).frame, 'r')
                fprintf("   Annotation at frame % 5d is: '%s'\n", annotations(i).frame, annotations(i).text{1})
            end
        end
    end
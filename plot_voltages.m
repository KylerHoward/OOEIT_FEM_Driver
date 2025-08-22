clearvars -except filepath
clc
close all

if exist('filepath', 'var') && ischar(filepath)
    [filename, filepath] = uigetfile(filepath, "Open Voltage File to Inspect");
else
    [filename, filepath] = uigetfile("Open Voltage File to Inspect");
end

fprintf("Inspecting %s\n", filename)
load(fullfile(filepath, filename))  
% Tripped at about 2190, so could only compute ECG to 2190
% Truncate the number of frames
frame_voltage = squeeze(frame_voltage(1:31,:,:));  % CP, electrode, frame
[tt,L,Slides] = size(frame_voltage);
fprintf("The voltage is %d by %d by %d\n", tt, L, Slides)

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
% plot them
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
    power_wavef = squeeze(power_wavef_mx(1,1,:));
    plot(power_wavef)
    title('Power Waveform First data set')
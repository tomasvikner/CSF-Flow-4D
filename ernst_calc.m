% Parameters
T1 = 1000;      % ms
TR = 12;        % ms  <-- change as needed

% Ernst angle (radians)
alpha_rad = acos(exp(-TR/T1));

% Convert to degrees
alpha_deg = alpha_rad * 180/pi;

fprintf('Ernst angle = %.2f degrees\n', alpha_deg);

%%
clc
close all

% Parameters
T1 = 1000;   % ms
TR = 10;     % ms

E1 = exp(-TR/T1);

% Flip angle range
FA = linspace(0.1,90,1000);    % degrees
FA_rad = FA*pi/180;

% SPGR signal (proportional to SNR)
S = (sin(FA_rad).*(1 - E1)) ./ (1 - cos(FA_rad).*E1);

% Normalize (optional)
S = S ./ max(S);

% Ernst angle
alphaE = acos(E1) * 180/pi;

% Plot
figure;
plot(FA, S, 'LineWidth',2);
hold on
xline(alphaE,'--');
xlabel('Flip Angle (deg)');
ylabel('Relative SNR');
title('SPGR SNR vs Flip Angle');
grid on

fprintf('Ernst angle = %.2f degrees\n', alphaE);
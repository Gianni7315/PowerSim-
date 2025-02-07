clear all;
close all;
clc;

% Design Specifications
fr = 100e3; % Resonant frequency in Hz
Vin_nom = 23; % Nominal input voltage in V
Vout = 22; % Output voltage in V
Nps = (Vin_nom / (2*Vout)); % Transformer turns ratio Half Bridge 
% for Full Bridge Nps = (Vin_nom / Vout);
Pout_max = 330; % Maximum output power in W
m = 6; % Ratio (Lr + Lm) / Lr
Q_max = .8; % Maximum quality factor

% Assume primary turns for practical implementation
Np = 1; % Example primary turns
Ns = round(Np / Nps); % Calculate secondary turns

% Calculate Reflected Load Resistance
Rac_min = (8 / (pi^2)) * (Nps^2 * Vout^2 / Pout_max);

% Calculate Resonant Tank Components
Lr = Rac_min * Q_max / (2 * pi * fr);
Cr = 1 / ((2 * pi * fr)^2 * Lr);
Lm = m * Lr - Lr;

% Calculate Output Load Resistance
Rout = Vout^2 / Pout_max;

% Normalized Frequency vs Gain Plot
F = linspace(0.1, 10, 10000); % Normalized frequency range
Q_values = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];

figure;
hold on;
for i = 1:length(Q_values)
    Q = Q_values(i);
    M = F.^2 * (m - 1) ./ sqrt((m * F.^2 - 1).^2 + F.^2 .* (F.^2 - 1).^2 * (m - 1)^2 * Q^2);
    plot(F, M, 'LineWidth', 1.5, 'DisplayName', sprintf('Q=%.2f', Q));
end

% Highlight the peak for Q_max
M_Qmax = F.^2 * (m - 1) ./ sqrt((m * F.^2 - 1).^2 + F.^2 .* (F.^2 - 1).^2 * (m - 1)^2 * Q_max^2);
[maxM, idx] = max(M_Qmax);
Fpeak = F(idx);

plot([Fpeak Fpeak], [0 maxM], 'r--', 'HandleVisibility', 'off');
plot(Fpeak, maxM, 'ro', 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');

xlabel('Normalized Frequency F_x');
ylabel('Gain M');
xlim([0.1 10]);
ylim([0 3]);
set(gca, 'XScale', 'log');
grid on;
legend show;
title('Normalized Frequency vs Gain (Multiple Q Values)');
hold off;

% Display computed values
disp('Resonant Tank Values:');
disp(['Lr = ', num2str(Lr*1e6), ' uH']);
disp(['Cr = ', num2str(Cr*1e6), ' uF']);
disp(['Lm = ', num2str(Lm*1e6), ' uH']);
disp(['Turns Ratio (Nps) = ', num2str(Nps)]);
disp(['Primary Turns (Np) = ', num2str(Np)]);
disp(['Secondary Turns (Ns) = ', num2str(Ns)]);
disp(['Output Load Resistance (Rout) = ', num2str(Rout), ' Ohms']);


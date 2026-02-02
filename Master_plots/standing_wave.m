% Standing Wave Visualization
clear; clc;
close all 

% Parameters
L = 1;              % Length of the string [m]
A = 1;              % Amplitude
v = 1;              % Wave speed [m/s]
n = 3;              % Mode number (1 = fundamental, 2 = first overtone, etc.)
T = 4;              % Total time of simulation [s]
fps = 120;           % Frames per second for animation

% Derived quantities
k = n * pi / L;                 % Wavenumber
omega = v * k;                  % Angular frequency
x = linspace(0, L, 500);        % Spatial grid
tVec = linspace(0, T, T*fps);   % Time grid

% Prepare figure
figure(1);
h = plot(x, zeros(size(x)), 'LineWidth', 2);
ylim([-1.2 1.2]*A);
xlabel('x [m]');
ylabel('Displacement');
title(['Standing Wave Mode n = ' num2str(n)]);
grid on;

% Animation loop
for t = tVec
    y = 2*A * sin(k*x) .* cos(omega*t); % Standing wave expression
    set(h, 'YData', y);
    drawnow;
end

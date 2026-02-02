% Parameters
Bo = 1/0.08;       % Bond number
Ma = 1;       % Marangoni number
Ga = 1;       % Galileo number
gamma0 = 2.7*1e-6;   % Reference shear

% Independent variable range
gamma = linspace(0, 5*gamma0, 500);

% Function definition
sigmaf = (1./Bo) .* (1 - tanh(Ma  ./ Ga .* (gamma ./ gamma0)));

% Plot
figure;
plot(gamma, sigmaf, 'b', 'LineWidth', 2);
xlabel('$\Gamma$','Interpreter','latex','FontSize',22);
ylabel('$\sigma_f$','Interpreter','latex','FontSize',22);
% title('\sigma_f vs. \gamma');
grid on;

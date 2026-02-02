 %% Master code to make all the plots that are going into the journal paper
qqq
%% Surface roughness plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%opening files

% x_temp = readmatrix('data/x_temp.txt');
% temp1= readmatrix('data/temp1.txt');
% xrange = readmatrix('data/xrange.csv');
% metadata.XResolution= readmatrix('data/XResolution.csv');
% yrange= readmatrix('data/yrange.csv');
% data_matrix= readmatrix('data/data_matrix.csv');
% y_temp= readmatrix('data/y_temp.txt');
% temp2= readmatrix('data/temp2.txt');
% xrange1= readmatrix('data/xrange1.csv');
% yrange1= readmatrix('data/yrange1.csv');
% wave_length1 = readmatrix('data/wave_length1.txt');
% wave_length2 = readmatrix('data/wave_length2.txt');
% %%Ploting
% figure(1);
% asetWH(gcf,10,9)
% subplot(2, 2, 1);
% plot(x_temp,temp1,'ko-','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','r');
% xlabel('$x \, (\mu m)$', 'Interpreter', 'latex', 'FontSize', 28);
% ylabel('$z \, (\mu m)$', 'Interpreter', 'latex', 'FontSize', 28);
% set(gca, 'FontSize', 28,'FontWeight','normal'); % Adjust the value as needed
% % Add annotation for dominant wavelength
% text_position_x = x_temp(round(length(x_temp)/2)-5); % Middle of the x_temp range
% text_position_y = max(temp1)-0.08; % Near the maximum of temp
% % text(text_position_x, text_position_y, ...
% %     sprintf('Dominant wavelength f: %.2f µm', wave_length1), ...
% %     'FontSize', 14, 'Color', 'r');
% axis padded
% asetWH(gcf,10,9)
% legend off
% % waterfall plot
% subplot(2, 2, 3);
% yrange = xrange;
% p = waterfall(xrange*metadata.XResolution,yrange*metadata.XResolution,data_matrix(xrange,yrange));
% p.EdgeColor = 'k';
% p.LineStyle = "-";
% p.Marker = '*';
% p.MarkerSize = 3;
% p.MarkerEdgeColor = 'r';
% xlabel('$x \, (\mu m)$', 'Interpreter', 'latex', 'FontSize', 28);
% ylabel('$y \, (\mu m)$', 'Interpreter', 'latex', 'FontSize', 28);
% zlabel('$z \, (\mu m)$', 'Interpreter', 'latex', 'FontSize', 28);
% asetWH(gcf,10,9)
% legend off
% subplot(2, 2, 2);
% plot(y_temp,temp2,'ko-','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','r');
% xlabel('$y \, (\mu m)$', 'Interpreter', 'latex', 'FontSize', 28);
% ylabel('$z \, (\mu m)$', 'Interpreter', 'latex', 'FontSize', 28);
% set(gca, 'FontSize', 28,'FontWeight','normal'); % Adjust the value as needed
% % Add annotation for dominant wavelength
% text_position_x = y_temp(round(length(y_temp)/2)-18); % Middle of the x_temp range
% text_position_y = max(temp2)-0.10; % Near the maximum of temp
% % text(text_position_x, text_position_y, ...
% % sprintf('Dominant wavelength : %.2f µm', wave_length2), ...
% % 'FontSize', 14, 'Color', 'r');
% axis padded
% asetWH(gcf,10,9)
% legend off
% % waterfall plot
% subplot(2, 2,4);
% xrange1 = yrange1;
% p = waterfall(xrange1*metadata.XResolution,yrange1*metadata.XResolution,data_matrix(xrange1,yrange1));
% set(gca, 'FontSize', 18,'FontWeight','normal'); % Adjust the value as needed
% axis normal
% p.EdgeColor = 'k';
% p.LineStyle = "-";
% p.Marker = '*';
% p.MarkerSize = 3;
% p.MarkerEdgeColor = 'r';
% % Add axis labels with TeX interpreter
% xlabel('$y \, (\mu m)$', 'Interpreter', 'latex', 'FontSize', 28);
% ylabel('$x \, (\mu m)$', 'Interpreter', 'latex', 'FontSize', 28);
% zlabel('$z \, (\mu m)$', 'Interpreter', 'latex', 'FontSize', 28);
% asetWH(gcf,10,9)
% legend off


%% Re V/S Ca plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the data
% Re V/S Ca plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read the data
Re_split = readmatrix('data/Re_split.csv');Ca_split = readmatrix('data/Ca_split.csv');
Re_rolling = readmatrix('data/Re_rolling.csv');Ca_rolling = readmatrix('data/Ca_rolling.csv');
% Define colors (RGB)
cl = [180, 210, 230;  % Pale Matte Blue (softer)
    230, 180, 180;  % Pale Matte Red (softer)
    200, 200, 200] / 255; % Pale Matte Grey (softer)
% Create figure
figure(2);
loglog(Re_split, Ca_split, 'pb', 'MarkerSize', 10, 'MarkerFaceColor', 'b'); hold on;
loglog(Re_rolling, Ca_rolling, 'or', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); hold on;
% Rectangle 1: Blue region at the bottom
x1 = [0.5, 0.5, 1e3, 1e3];
y1 = [1e-6, 2e-5, 2e-5, 1e-6];
fill(x1, y1, cl(1,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
% Rectangle 2: Gradient region (Blue → Red from bottom to top)
x2 = [0.5, 1e3, 1e3, 0.5];
y2 = [2e-5, 2e-5, 2e-4, 2e-4];
% Define colors at each vertex to achieve vertical gradient
cmap = [cl(1,:); cl(1,:); cl(2,:); cl(2,:)]; % Bottom two blue, top two red
% Use `patch` for smooth interpolation from bottom to top
p = patch(x2, y2, zeros(size(x2)), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
set(p, 'FaceVertexCData', cmap, 'FaceColor', 'interp');
% Rectangle 3: Red region at the top
x3 = [0.5, 1e3, 1e3, 0.5];
y3 = [2e-4, 2e-4, 1, 1];
fill(x3, y3, cl(2,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
% Grid, labels, and legend
grid on;
xlabel('$Re_f$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Ca', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold', 'Rotation', 90);
legend({'Split Streamline', 'Rolling'}, 'Location', 'best');
axis tight
% Adjust figure size
asetWH(gcf, 6, 5);
% % %%%% Construct the new filename
% print('The plot is to be saved as a pdf and .fig file')
% fname = 'ReVSCa.pdf';
% figfname = 'ReVSCa.fig';
% savefig(figfname);
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print(gcf, '-dpdf','-r600', fname);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Interface shape comparison plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading interface data
f8 = readmatrix("data/f_avg8.txt");f9 = readmatrix("data/f_avg9.txt");
f10 = readmatrix("data/f_avg10.txt"); f11 = readmatrix("data/f_avg11.txt");
x11=readmatrix('data/x11.txt');y11=readmatrix('data/y11.txt');
x10=readmatrix('data/x10.txt'); y10=readmatrix('data/y10.txt');
x9=readmatrix('data/x9.txt');y9=readmatrix('data/y9.txt');
x8=readmatrix('data/x8.txt');y8=readmatrix('data/y8.txt');
% Define colors as RGB triplets
cm = [ 1/255, 31/255, 75/255;    % Dark Navy Blue (#011f4b)
    3/255, 57/255, 108/255;   % Deep Blue (#03396c)
    0/255, 91/255, 150/255;   % Strong Blue (#005b96)
    100/255, 151/255, 177/255; % Soft Blue (#6497b1)
   ];
% To find the location of the interface
idx10 = find(0.45<f10 & 0.55>f10 );
idx11 = find(0.45<f11 & 0.55>f11 );
idx9 = find(0.45<f9 & 0.55>f9 );
idx8 = find(0.2<f8 & 0.8>f8 );
figure(3)
plot(x11(idx11), y11(idx11), '--*', 'Color', cm(4, :), 'LineWidth', 2);
hold on
plot(x10(idx10), y10(idx10), '--*', 'Color', cm(3, :), 'LineWidth', 2);
plot(x9(idx9), y9(idx9), '--*', 'Color', cm(2, :), 'LineWidth', 2);
plot(x8(idx8), y8(idx8), '--*', 'Color', cm(1, :), 'LineWidth', 2);
asetWH(gcf,6,6)
xlabel('x(m)');
ylabel('y(m)');
legend('11','10','9','8')
hold on
grid on
asetWH(gcf,6,6)
% axis equal
% %Construct the new filename
% [~, current_folder_name] = fileparts(pwd);
% fname = sprintf('interface_shape_comparison%s.pdf', current_folder_name);
% figfname = sprintf('interface_shape_comparison%s.fig', current_folder_name);
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print(gcf, '-dpdf','-r600', fname);
% % Save the file
% savefig(figfname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% grid independence study
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the data from files
data5 = readtable("data/data5.txt");data11 = readtable("data/data11.txt");
data7 = readtable("data/data7.txt");data6 = readtable("data/data6.txt");
data9 = readtable("data/data9.txt");data8 = readtable("data/data8.txt");
data10 = readtable("data/data10.txt");
U0 = 0.0001;L0 = 0.025;f = U0 / L0;
% Extract x and y values for each dataset
x11 = f*table2array(data11(:,"Var2")); y11 = table2array(data11(:,"Var6"));
x5 = f*table2array(data5(1:length(x11),"Var2")); y5 = table2array(data5(1:length(y11),"Var6"));
x6 = f*table2array(data6(1:length(x11),"Var2")); y6 = table2array(data6(1:length(y11),"Var6"));
x7 = f*table2array(data7(1:length(x11),"Var2")); y7 = table2array(data7(1:length(y11),"Var6"));
x8 = f*table2array(data8(1:length(x11),"Var2")); y8 = table2array(data8(1:length(y11),"Var6"));
x9 = f*table2array(data9(1:length(x11),"Var2")); y9 = table2array(data9(1:length(y11),"Var6"));
x10 = f*table2array(data10(1:length(x11),"Var2")); y10 = table2array(data10(1:length(y11),"Var6"));
% Common x-values for calculating differences
x_common = x10; % Assuming all datasets share the same x-values range.
% Main Plot
lw = 'LineWidth'; ms = 'MarkerSize';mfc = 'MarkerFaceColor';mec = 'MarkerEdgeColor';
line_colors = [
    0/255, 76/255, 76/255;     % Deep Teal (#004c4c)
    0/255, 102/255, 102/255;   % Dark Teal (#006666)
    0/255, 128/255, 128/255;   % Teal (#008080)
    102/255, 178/255, 178/255; % Soft Teal-Green (#66b2b2)
    150/255, 200/255, 200/255; % Light Cyan-Green (#b2d8d8)
    178/255, 216/255, 216/255; % Light Cyan-Green (#b2d8d8)
];
marker_types = {'-d','-v',  '-p', '-^ ' , '-h', '-o', '-s'};
prange = 1:1:1100;
prange1 = 1:450:1100;
figure(4);
hold on;
% Plot the main data
plot(x6(prange), y6(prange), lw, 1.7);
plot(x7(prange), y7(prange),  lw, 1.7);
plot(x8(prange), y8(prange), lw, 1.7);
plot(x9(prange), y9(prange), lw, 1.7);
plot(x10(prange), y10(prange), lw, 1.7);
plot(x11(prange), y11(prange), lw, 1.7);
colororder(line_colors)
xlim([0 0.022])
ylim([0 0.0022])
hold on
%
plot(x6(1:390:1100), y6(1:390:1100), 'LineStyle', 'none', 'Marker', marker_types{7}(2:end), ms, 12, mec, 'k', mfc, line_colors(1,:));
plot(x7(1:400:1100), y7(1:400:1100), 'LineStyle', 'none', 'Marker', marker_types{5}(2:end), ms, 12, mec, 'k', mfc, line_colors(2,:));
plot(x8(1:420:1100), y8(1:420:1100), 'LineStyle', 'none', 'Marker', marker_types{4}(2:end), ms, 12, mec, 'k', mfc,  line_colors(3,:));
plot(x9(1:450:1100), y9(1:430:1100), 'LineStyle', 'none', 'Marker', marker_types{3}(2:end), ms, 12, mec, 'k', mfc,  line_colors(4,:));
plot(x10(1:470:1100), y10(1:440:1100), 'LineStyle', 'none', 'Marker', marker_types{2}(2:end), ms, 12, mec, 'k', mfc, line_colors(5,:));
plot(x11(1:490:1100), y11(1:490:1100), 'LineStyle', 'none', 'Marker', marker_types{1}(2:end), ms, 12, mec, 'k', mfc,  line_colors(6,:));
hold off;
grid on;
% Labels and legend
xlabel('$\tau$', 'Interpreter', 'latex');
ylabel('$h_{CL} (m)$', 'Interpreter', 'latex');
legend({ '$\Delta = 1/64$', '$\Delta = 1/128$', '$\Delta = 1/256$', ...
    '$\Delta = 1/512$', '$\Delta = 1/1024$', '$\Delta = 1/2048$'}, ...
    'Location', 'best', 'Interpreter', 'latex');
asetWH(gcf, 6, 5);
legend('FontSize', 16, 'NumColumns', 1,'Location','southwest');
% Define inset position (normalized units: [left bottom width height])
inset_position = [0.57, 0.25, 0.3, 0.3];
% Create inset axes
ax_inset = axes('Position', inset_position);
hold on;
Delta = [1/2^11, 1/2^10, 1/2^9, 1/2^8, 1/2^7, 1/2^6 ];
hcl_end_steady = [y11(end), y10(end), y9(end), y8(end), y7(end), y6(end) ];
% Plot the inset
semilogx(Delta, hcl_end_steady,'-o', 'Color', [0/255, 76/255, 76/255], 'MarkerSize', 10, 'LineWidth', 1.5, mfc, 'w');
grid off;	
axis padded;
% Customize inset appearance
set(ax_inset, 'XColor', 'k', 'YColor', 'k', 'Box', 'on', 'FontSize', 10);
xlabel(ax_inset, '$\Delta$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel(ax_inset, '$h_{CL} (m)$', 'Interpreter', 'latex', 'FontSize', 12);
hold off;
% %%%% Construct the new filename
print('The plot is to be saved as a pdf and .fig file')
fname = 'gridindependence.pdf';
figfname = 'gridindependence.fig';
savefig(figfname);
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print(gcf, '-dpdf','-r600', fname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Load the saved plotting data
load('./data/plotting_ROLLING_data.mat');

% Apply the lowpass filter again if not already saved
% (Uncomment this if 'real_psi' wasn't saved)
% real_psi = lowpass_fil(psi, 300);
figure;
stp = stp*2;
level = 1*[0.01:0.05:0.9];
label_levels = [ 0 0.2 0.6 ];
hold on
u_mag = sqrt(u_m.^2 + v_m.^2);
h1 = pcolor(x_m, y_m, u_mag); 
shading interp;
colormap(colorcet('COOLWARM'));
colorbar
clim([0 max(u_mag(:))]);
[M, c] = contour(x_m, y_m, real_psi, 'LineWidth', 2.5);
% clabel(M, c,  label_levels, 'FontWeight', 'bold', 'Color', 'w', ...
       % 'FontSize', 12, 'Interpreter', 'latex', 'Backgroundcolor', 'none');
quiver(x_m(1:stp:end, 1:stp:end), y_m(1:stp:end, 1:stp:end), ...
       u_m(1:stp:end, 1:stp:end), v_m(1:stp:end, 1:stp:end), 'k', 'LineWidth', 1.2);
plot(x(idx), y(idx), 'k', 'LineWidth', 2.5);
xlabel('$x/l_c$', 'Interpreter', 'latex','fontsize',20);
ylabel('$y/l_c$', 'Interpreter', 'latex','fontsize',20);

axis equal tight
xlim([0 max(x_m(:))/3.5]);
ylim([min(y_m(:))/3 max(y_m(:))/2]);
xticks(0:2e-3:15e-3);
yticks(-8e-3:2e-3:8e-3);
xticklabels(arrayfun(@(x) num2str(round(x/L)), xticks, 'UniformOutput', false));
yticklabels(arrayfun(@(y) num2str(round(y/L)), yticks, 'UniformOutput', false));
set(gca, 'TickLength', [0.02, 0.0], 'XMinorTick', 'off', 'YMinorTick', 'off');
set(gca, 'FontSize', 16,'TickLabelInterpreter', 'latex')
hold off



%%
clear 
load('./data/plotting_SPLIT_data.mat');
figure(2)
stp = stp*2;
level = 1*[0.01:0.05:0.9];
label_levels = [ 0 0.2 0.6 ];
hold on
u_mag = sqrt(u_m.^2 + v_m.^2);
h1 = pcolor(x_m, y_m, u_mag); 
shading interp;
colormap(colorcet('COOLWARM'));
colorbar
clim([0 max(u_mag(:))]);
[M, c] = contour(x_m, y_m, real_psi, 'LineWidth', 2.5);
% clabel(M, c,  label_levels, 'FontWeight', 'bold', 'Color', 'w', ...
       % 'FontSize', 12, 'Interpreter', 'latex', 'Backgroundcolor', 'none');
quiver(x_m(1:stp:end, 1:stp:end), y_m(1:stp:end, 1:stp:end), ...
       u_m(1:stp:end, 1:stp:end), v_m(1:stp:end, 1:stp:end), 'k', 'LineWidth', 1.2);
plot(x(idx), y(idx), 'k', 'LineWidth', 2.5);
xlabel('$x/l_c$', 'Interpreter', 'latex','fontsize',20);
ylabel('$y/l_c$', 'Interpreter', 'latex','fontsize',20);

axis equal tight
xlim([0 max(x_m(:))/3.5]);
ylim([min(y_m(:))/3 max(y_m(:))/2]);
xticks(0:2e-3:15e-3);
yticks(-8e-3:2e-3:8e-3);
xticklabels(arrayfun(@(x) num2str(round(x/L)), xticks, 'UniformOutput', false));
yticklabels(arrayfun(@(y) num2str(round(y/L)), yticks, 'UniformOutput', false));
set(gca, 'TickLength', [0.02, 0.0], 'XMinorTick', 'off', 'YMinorTick', 'off');
set(gca, 'FontSize', 16,'TickLabelInterpreter', 'latex')
% % Construct the new filename
fname = sprintf('./figs_final/split.pdf');
fname1 = sprintf('./figs_final/split.png');
figfname = sprintf('./figs_final/split.fig');
%  % Construct the new filename
savefig(figfname);
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
exportgraphics(gcf,fname,Resolution=600);
exportgraphics(gcf,fname1,Resolution=600);



%% %% flow fields gly Re9
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_avg= readmatrix('data/re9/f_avg.txt');psi= readmatrix('data/re9/psi.txt');
% X= readmatrix('data/re9/X.txt');Y= readmatrix('data/re9/Y.txt');
% v_bar= readmatrix('data/re9/v_bar.txt');u_bar= readmatrix('data/re9/u_bar.txt');
% x = readmatrix('data/re9/x.txt');y = readmatrix('data/re9/y.txt');
% stp = 20;
% % f_avg = mean(f_all,2);
% idx = find(0.1<f_avg & 0.9>f_avg );
% real_psi = lowpass_fil(psi,320);
% figure(1);
% % levels = [-logspace(-4, 0, 10) logspace(-4, 0, 10)];
% levels = [linspace(-0.3,0,5) linspace(0,1.5,10)];
% hold on
% contourf(X, Y, real_psi, levels,'W', 'LineWidth', 1.5);
% hold on
% plot(x(idx),y(idx),'r',LineWidth=3.5)
% quiver(X(1:stp:end, 1:stp:end), Y(1:stp:end, 1:stp:end), u_bar(1:stp:end, 1:stp:end), ...
%     v_bar(1:stp:end, 1:stp:end),'k','LineWidth',1.2);
% % plot(x_filtered, y_filtered, 'r', 'LineWidth', 2.5);  % Interface positions
% colorbar
% colorcet('COOLWARM')
% hold off
% axis equal;
% axis tight;
% % Get current tick values
% xticks(0:2e-3:15*1e-3)
% yticks(-8*1e-3:2e-3:8*1e-3)
% xt = xticks; % Get x tick values
% yt = yticks; % Get y tick values
% % Convert meters to millimeters
% xticklabels(num2str(xt(:) * 1000)); % Convert to mm and format
% yticklabels(num2str(yt(:) * 1000)); % Convert to mm and format
% % Add units in the axis labels
% xlabel('x (mm)')
% ylabel('y (mm)')
% asetWH(gcf,6,6)
% set(gca, 'TickLength', [0.02, 0.0]); % Major tick: 0.02, Minor tick: 0.01
% set(gca, 'XMinorTick', 'off', 'YMinorTick', 'off'); % Turn on minor ticks for both axes
% legend off
% xlim([0 max(max(X))/3]);
% ylim([min(min(Y))/3 max(max(Y))/2.5]);
% % Construct the new filename
% fname = sprintf('sl_fftfilter_closeup_re9.pdf');
% figfname = sprintf('sl_fftfilter_closeup_re9.fig');
% %  % Construct the new filename
% print('The plot is to be saved as a pdf and .fig file')
% savefig(figfname);
% 
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% exportgraphics(gcf,fname,Resolution=600);
% %
% figure(2);
% % levels = [-logspace(-4, 0, 10) logspace(-4, 0, 10)];
% % levels = [linspace(-0.3,0,5) linspace(0,1.5,10)];
% levels = [ logspace(-5,-3.5,4) logspace(-2.5,-1,4) linspace(0.1,1,9) -linspace(0.1,0.25,7) -0.045 -0.016 ...
%     -0.0019 -0.0058];
% r_psi=lowpass_fil(psi,50);
% hold on
% contourf(X, Y, r_psi, levels,'w', 'LineWidth', 0.8);
% hold on
% plot(x(idx),y(idx),'r',LineWidth=1.5)
% % quiver(X(1:stp:end, 1:stp:end), Y(1:stp:end, 1:stp:end), u_bar(1:stp:end, 1:stp:end), ...
% %     v_bar(1:stp:end, 1:stp:end),'m','LineWidth',1.2);
% % plot(x_filtered, y_filtered, 'r', 'LineWidth', 2.5);  % Interface positions
% colorbar
% colorcet('COOLWARM')
% hold off
% axis equal;
% axis tight;
% % Get current tick values
% xticks(0:2e-3:15*1e-3)
% yticks(-8*1e-3:2e-3:8*1e-3)
% xt = xticks; % Get x tick values
% yt = yticks; % Get y tick values
% % Convert meters to millimeters
% xticklabels(num2str(xt(:) * 1000)); % Convert to mm and format
% yticklabels(num2str(yt(:) * 1000)); % Convert to mm and format
% % Add units in the axis labels
% xlabel('x (mm)')
% ylabel('y (mm)')
% asetWH(gcf,6,6)
% ax = gca;
% ax.YAxis.Exponent = -3;
% ax.XAxis.Exponent = -3;
% set(gca, 'TickLength', [0.02, 0.0]); % Major tick: 0.02, Minor tick: 0.01
% set(gca, 'XMinorTick', 'off', 'YMinorTick', 'off'); % Turn on minor ticks for both axes
% legend off
% 
% % Construct the new filename
% fname = sprintf('sl_fftfilter_re9.pdf');
% figfname = sprintf('sl_fftfilter_re9.fig');
% %  % Construct the new filename
% print('The plot is to be saved as a pdf and .fig file')
% savefig(figfname);
% 
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% exportgraphics(gcf,fname,Resolution=600);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% New sims Re V/S Ca plot 
clc; clear; close all

% Split data
% ReF_split = [3.021164;14.699598;18.374498;18.374498;20.871402;83.326763;29.399197;183.744981;
%     33.412051;186.890761;186.890761;208.714017;711.401214;833.267630;1043.570087;2087.140173];
% Ca_split = [1.39E-05;1.39E-05;6.94E-06;1.39E-05;2.78E-06;2.78E-06;1.39E-05;6.94E-06;
%     1.39E-05;2.78E-06;1.39E-05;1.39E-05;1.39E-06;1.39E-05;1.39E-05;1.39E-05];
ReF_split = [147.556891;14.699598;18.374498;18.374498;20.871402;83.326763;29.399197;183.744981;
    33.412051;186.890761;186.890761;208.714017;711.401214;833.267630;1043.570087;2087.140173];
Ca_split = [9.82E-04;1.39E-05;6.94E-06;1.39E-05;2.78E-06;2.78E-06;1.39E-05;6.94E-06;
    1.39E-05;2.78E-06;1.39E-05;1.39E-05;1.39E-06;1.39E-05;1.39E-05;1.39E-05];

% Rolling data
ReF_rolling = [9.048256;0.367490;0.734980;1.837450;1.837450;3.674900;1.812698;0.642398;25.647041;
    32.119914;64.239828;6.734468;13.109658;52.338858];
Ca_rolling = [1.60E-03;1.39E-05;1.39E-05;1.39E-03;1.39E-05;1.39E-03;1.39E-05;6.77E-03;4.51E-04;
    2.26E-03;2.26E-03;6.72E-03;4.42E-03;2.21E-01];
pt_rolReF = 9.048256;
pt_rolCa= 1.60E-03;
pt_splitReF = 147.556891;
pt_splitCa= 9.82E-04;
% Combine data
X = log10([ReF_split; ReF_rolling]); % log ReF
Y = log10([Ca_split; Ca_rolling]);   % log Ca
data = [X, Y];
labels = [ones(length(ReF_split), 1); -ones(length(ReF_rolling), 1)]; % 1: split, -1: rolling

% Train linear SVM
SVMModel = fitcsvm(data, labels, 'KernelFunction','polynomial');

% Get decision boundary
w = SVMModel.Beta; b = SVMModel.Bias;

% Mesh for decision region visualization
[xGrid, yGrid] = meshgrid(linspace(min(X)-0.1, max(X)+0.1, 500), ...
                          linspace(min(Y)-1, max(Y)+1, 500));
gridPoints = [xGrid(:), yGrid(:)];
[~, score] = predict(SVMModel, gridPoints);
scoreGrid = reshape(score(:,2), size(xGrid));

% Plot background regions
figure; hold on;
contourf(10.^xGrid, 10.^yGrid, scoreGrid > 0, 1, 'LineColor', 'none');
colormap([1 0.85 0.85; 0.85 0.9 1]); % Rolling = redish, Split = blueish

% Plot decision boundary
% Plot the thick gray line
[C, hContour] = contour(10.^xGrid, 10.^yGrid, scoreGrid, [0 0], ...
    'Color', [0.9 0.9 0.9], 'LineWidth', 12);  % gray color
[C, hContour] = contour(10.^xGrid, 10.^yGrid, scoreGrid, [0 0], ...
    'Color', [0.8 0.8 0.8], 'LineWidth', 11);  % gray color
[C, hContour] = contour(10.^xGrid, 10.^yGrid, scoreGrid, [0 0], ...
    'Color', [0.7 0.7 0.7], 'LineWidth', 10);  % gray color
[C, hContour] = contour(10.^xGrid, 10.^yGrid, scoreGrid, [0 0], ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 9);  % gray color




% Plot actual data points
markerSize = 16;

% Plot actual data points using plot
hold on
s1 = plot(ReF_split, Ca_split, 'p', ...
    'MarkerSize', markerSize, ...
    'MarkerEdgeColor', 'b',MarkerFaceColor='b'); % Star for Split

s2 = plot(ReF_rolling, Ca_rolling, 'o', ...
    'MarkerSize', markerSize, ...
    'MarkerEdgeColor', 'r', ...
    'MarkerFaceColor', 'r'); % Filled circle for Rolling
 plot(pt_rolReF,pt_rolCa,'ok','MarkerSize', 10,'MarkerEdgeColor', 'k',MarkerFaceColor='w')
 plot(pt_splitReF,pt_splitCa,'pk','MarkerSize', 10,'MarkerEdgeColor', 'k',MarkerFaceColor='w')
 hold off


% Axes formatting
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('$Re_f$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('$Ca$', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
legend( '','','','','','Rolling away ','Rolling in', ...
       'Location','northeast');
% title('Flow Regime Demarcation by SVM')
asetWH(gcf,7,6); % Optional size control




%% Flow fields for rolling and split pattern 
close all; clear; clc

figure(5); clf;
fs=24;

% ================= ROLLING =================
load('./data/plotting_ROLLING_data_RF11.mat');
figure(5);
level = 1*[0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,...
    0.65,0.7,0.75,0.8,0.85,0.9];
hold on
u_mag = sqrt(u_m.^2+v_m.^2);
pcolor(x_m/lc,y_m/lc,u_mag), hold on;colorcet('COOLWARM'); shading interp;
colorbar, clim([0 max(max(u_mag))]);

[M,c] = contour(x_m/lc, y_m/lc, psi,'k', 'LineWidth', 2.5);

quiver(x_m(1:stp:end, 1:stp:end)/lc, y_m(1:stp:end, 1:stp:end)/lc, u_m(1:stp:end, 1:stp:end), ...
    v_m(1:stp:end, 1:stp:end),'k','LineWidth',0.6);

plot(x(idx)/lc, y(idx)/lc, 'b', 'LineWidth', 3.5); % Interface positions
axis equal;
axis tight;
hold off
% Add units in the axis labels
xlabel('$ x/l_c$ ','Interpreter','latex','FontSize',fs)
ylabel('$ y/l_c$ ','Interpreter','latex','FontSize',fs)
set(gca,'fontsize',fs,'FontName','latex')
% Get the current axes handle
ax = gca;
% Set the interpreter for the tick labels to LaTeX
ax.TickLabelInterpreter = 'latex';
% NOTE: The asetWH function is not a standard MATLAB function.
% You will need to define it or remove it.
% asetWH(gcf,6,6)
set(gca, 'TickLength', [0.02, 0.0]); % Major tick: 0.02, Minor tick: 0.01
set(gca, 'XMinorTick', 'off', 'YMinorTick', 'off'); % Turn on minor ticks for both axes
xlim([0 max(max(x_m/lc))/3]);
ylim([min(min(y_m/lc))/2.5 max(max(y_m/lc))/2.]);
grid off
asetWH(gcf, 6,7);legend off

% ================= SPLIT =================
load('./data/plotting_SPLIT_data_RF11.mat');
figure(6);asetWH(gcf, 7, 6);
level = 1*[0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,...
    0.65,0.7,0.75,0.8,0.85,0.9];
hold on
u_mag = sqrt(u_m.^2+v_m.^2);
bgcolor = u_mag;
pcolor(x_m/lc,y_m/lc,bgcolor), hold on;colorcet('COOLWARM'); shading interp;
% colorbar, clim([(mean(mean(bgcolor))- 1*std(std(bgcolor)))  (mean(mean(bgcolor))+ 1*std(std(bgcolor))) ]);
colorbar, clim([0 max(max(bgcolor)) ]);

xlim([0 max(max(x_m/lc))/3]);
ylim([min(min(y_m/lc))/2.5 max(max(y_m/lc))/2.]);

[M,c] = contour(x_m/lc, y_m/lc, psi,'k', 'LineWidth', 2.5);

quiver(x_m(1:stp:end, 1:stp:end)/lc, y_m(1:stp:end, 1:stp:end)/lc, u_m(1:stp:end, 1:stp:end), ...
    v_m(1:stp:end, 1:stp:end),'k','LineWidth',0.6);

plot(x(idx)/lc, y(idx)/lc, 'b', 'LineWidth', 3.5); % Interface positions
% Extract the interface points
x_int = x(idx);
y_int = y(idx);

axis equal;
axis tight;
hold off
% Add units in the axis labels
xlabel('$ x/l_c$ ','Interpreter','latex','FontSize',fs)
ylabel('$ y/l_c$ ','Interpreter','latex','FontSize',fs)
set(gca,'fontsize',fs,'FontName','latex')
% Get the current axes handle
ax = gca;
% Set the interpreter for the tick labels to LaTeX
ax.TickLabelInterpreter = 'latex';
% NOTE: The asetWH function is not a standard MATLAB function.
% You will need to define it or remove it.
% asetWH(gcf,6,6)
set(gca, 'TickLength', [0.02, 0.0]); % Major tick: 0.02, Minor tick: 0.01
set(gca, 'XMinorTick', 'off', 'YMinorTick', 'off'); % Turn on minor ticks for both axes
xlim([0 max(max(x_m/lc))/3]);
ylim([min(min(y_m/lc))/2.5 max(max(y_m/lc))/2.]);
asetWH(gcf, 6,7);legend off


%% This part of the code plots the Hysteresis loop, CL ANGLE,POSITION , VELOCITY evolution 


data = readmatrix("data/CL_info_rolling.txt");

f = 100; stp = 100; TP = 1/f;
t0 = data(:,2); t1 = data(:,2) - data(1,2); t = t1 / TP;
ycl = data(:,4); Ucl_mag = data(:,5); Ucl_u = data(:,6); Ucl_v = data(:,7);
theta_var = rad2deg(data(:,8));

% Normalized quantities
NORM_ycl = (ycl - mean(ycl)) / max(ycl - mean(ycl));
U_mag = sign(Ucl_u) .* sqrt(Ucl_u.^2 + Ucl_v.^2);
NDU_mag = (U_mag - mean(U_mag)) / max(U_mag - mean(U_mag));
ND_theta = theta_var / abs(max(theta_var));
NORM_theta = (theta_var - mean(theta_var)) / max(theta_var - mean(theta_var));
NORM_v = (Ucl_v - mean(Ucl_v)) / max(Ucl_v - mean(Ucl_v));

% Stroke and umax
format long
stroke = (max(ycl) - min(ycl)) * 1000; % in mm
umax = f * stroke / 1000;

% Plot: Phase-normalized signals
figure(1); hold on
fin = 200;

% --- Left Y-axis ---
yyaxis left
plot(t(1:fin), NORM_theta(1:fin), 'b-', 'LineWidth', 2);
plot(t(1:6:fin), NORM_theta(1:6:fin), 'bo', 'MarkerEdgeColor', 'k', ...
     'MarkerFaceColor', 'b', 'MarkerSize', 10);
plot(t(1:fin), NORM_ycl(1:fin), 'b-', 'LineWidth', 2);
plot(t(1:6:fin), NORM_ycl(1:6:fin), 'gv', 'MarkerEdgeColor', 'k', ...
     'MarkerFaceColor', 'g', 'MarkerSize', 10);
ylabel('$\phi$ \& $\hat{y_{\mathrm{cl}}}$', 'Interpreter', 'latex');
ax = gca; ax.YColor = 'b';

xlabel('t/T');

% --- Right Y-axis ---
yyaxis right
plot(t(1:fin), NORM_v(1:fin) / max(NORM_v), 'r-', 'LineWidth', 2);
plot(t(1:6:fin), NORM_v(1:6:fin) / max(NORM_v), 'rs', ...
     'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
% plot(t(1:100), NDU_mag(1:100), 'r-', 'LineWidth', 2);
% plot(t(1:6:100), NDU_mag(1:6:100), 'rs', ...
%      'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
ylabel('$\hat{U_{\mathrm{cl}}}$', 'Interpreter', 'latex');
ax.YAxis(2).Color = 'r';

% --- Styling ---
asetWH(gcf, 6, 3);
legend({'','$\phi$','','$\hat{y_{\mathrm{cl}}}$','', '$\hat{U_{\mathrm{cl}}}$'}, ...
       'Interpreter', 'latex', 'FontSize', 16);
grid on; hold off

% Save Figure 1
figure(1)
fname = './fig_int/HU_.png'; fname2 = './fig_int/HU_.pdf'; figfname = './fig_int/HU_.fig';
savefig(figfname);
set(gcf, 'Units', 'inches');
screenposition = get(gcf, 'Position');
set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], ...
         'PaperSize', screenposition(3:4));
exportgraphics(gcf, fname, 'Resolution', 600);
exportgraphics(gcf, fname2, 'Resolution', 600);

% Hysteresis loop
figure(2)
plot(-Ucl_v(1:275), theta_var(1:275), 'o-', ...
     'LineWidth', 1.5, 'MarkerFaceColor', 'w');
ylabel('$\phi$', 'Interpreter', 'latex');
xlabel('$U_{\mathrm{cl}}$', 'Interpreter', 'latex');
axis padded; asetWH(gcf, 4, 4); legend off

% Save Figure 2
% figure(2)
% fname = './fig_int/hysteresis_loop.png'; 
% fname2 = './fig_int/hysteresis_loop.pdf'; 
% figfname = './fig_int/hysteresis_loop.fig';
% savefig(figfname);
% set(gcf, 'Units', 'inches');
% screenposition = get(gcf, 'Position');
% set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], ...
%          'PaperSize', screenposition(3:4));
% exportgraphics(gcf, fname, 'Resolution', 600);
% exportgraphics(gcf, fname2, 'Resolution', 600);

% 3D Hysteresis Visualization
endin = 100;
displacement = NORM_ycl(1:endin);
velocity = Ucl_v(1:endin);
angle = theta_var(1:endin);

% % Optional smoother version (use instead of spline)
% windowSize = 9; % choose an odd number
% displacement = movmean(displacement, windowSize);
% velocity = movmean(velocity, windowSize);
% angle = movmean(angle, windowSize);

figure(3); hold on
% Index markers (if needed)
one = 1:25; two = 25:50; three = 50:75; four = 75:100;

% Main 3D curve
plot3(velocity, displacement, angle, 'o-', ...
       'MarkerSize', 5, 'MarkerFaceColor', [0.3 0.5 0.8], ...
       'Color', [0.2 0.4 0.6], 'LineWidth', 1.2);

xlabel('Contact line velocity (m s^{-1})');
ylabel('Contact line displacement (mm)');
zlabel('Contact angle (deg.)');
box on; axis padded; view(135, 50);

% --- Projections ---
plot3(velocity, displacement, 0.8*min(angle)*ones(size(angle)), '--', ...
      'Color', [0.6 0.6 0.1]);  % XY plane
plot3(velocity, 1.2*min(displacement)*ones(size(displacement)), angle, '--', ...
      'Color', [0.1 0.6 0.6]);  % XZ plane
plot3(1.2*min(velocity)*ones(size(velocity)), displacement, angle, '--', ...
      'Color', [0.6 0.1 0.6]);  % YZ plane

hold off;

% % Save 3D figure
% figure(1)
% fname = './fig_int/3D_plot_yuth.png';
% fname2 = './fig_int/3D_plot_yuth.pdf';
% figfname = './fig_int/3D_plot_yuth.fig';
% savefig(figfname);
% set(gcf, 'Units', 'inches');
% screenposition = get(gcf, 'Position');
% set(gcf, 'PaperPosition', [0 0 screenposition(3:4)], ...
%          'PaperSize', screenposition(3:4));
% exportgraphics(gcf, fname, 'Resolution', 600);
% exportgraphics(gcf, fname2, 'Resolution', 600);


%% Stationary plate contactline time evolving graphs
ms =12;lw=1.5;
%------------------------------------------------------------
%+++++++++++++++++++++++dt-CA90~~~~f=5
%~------------------------------------------------------------
fpath= "./data/CL_info_statPL_f5.txt";
data = readmatrix(fpath,'NumHeaderLines',1); 
% ------------------------------------------------------------
%  Data Processing
% ------------------------------------------------------------
f = 5;h=0.001345;phi_not=deg2rad(15);
TP = 1/f;t0 = data(:, 3); t1 = t0 - t0(1);
t  = t1 ./ TP;
ycl     = data(:, 5) ; 
Ucl_mag = data(:, 6);Ucl_u   = data(:, 7) ;
Ucl_v   = data(:, 8) ;theta_var = data(:, 9);
% ------------------------------------------------------------
% Normalisation (same logic as Python)
% ------------------------------------------------------------
NDU_v = (Ucl_v - mean(Ucl_v))/(h*f);
ND_ycl = (ycl - mean(ycl))/h;
NORM_theta = (theta_var - mean(theta_var))/phi_not;
theta_eq = mean(theta_var);
% % ------------------------------------------------------------
% Plot: Phase-normalized signals
% ------------------------------------------------------------
figure(1);asetWH(gcf, 8, 8);legend off; hold on
fin = 199;
% % --- Left Y-axis ---
subplot(3,1,1);hold on
plot(t(1:fin), NORM_theta(1:fin), 'w--', 'LineWidth', lw);
plot(t(1:5:fin), NORM_theta(1:5:fin), 'bs', 'MarkerEdgeColor', 'k', ...
     'MarkerFaceColor', 'k', 'MarkerSize', ms);
ylabel('$\phi$', 'Interpreter', 'latex');
axis tight
xlabel('t/T');
asetWH(gcf, 6, 3);legend off
grid on; hold off
subplot(3,1,2);hold on
plot(t(1:fin), ND_ycl(1:fin), 'w--', 'LineWidth', lw);
plot(t(1:6:fin), ND_ycl(1:6:fin), 'bs', 'MarkerEdgeColor', 'b', ...
     'MarkerFaceColor', 'b', 'MarkerSize', ms);
ylabel('$\hat{y_{\mathrm{cl}}}$', 'Interpreter', 'latex');
axis tight
xlabel('t/T');
asetWH(gcf, 6, 6);legend off
grid on; hold off
subplot(3,1,3);hold on
plot(t(1:fin), NDU_v(1:fin), 'w--', 'LineWidth', lw);
plot(t(1:6:fin), NDU_v(1:6:fin), 'rs', 'MarkerEdgeColor', 'r', ...
     'MarkerFaceColor', 'r', 'MarkerSize', ms);
ylabel('$\hat{U_{\mathrm{cl}}}$', 'Interpreter', 'latex');
axis tight
xlabel('t/T');
asetWH(gcf, 6, 3);legend off
grid on; hold off
%+++++++++++++++++++++++dt-15
fpath= "./data/CL_info_statPL_f100.txt";
data = readmatrix(fpath,'NumHeaderLines',1); 
% % % ------------------------------------------------------------
 % % % Data Processing
% % % % ------------------------------------------------------------
f = 100;h=0.001345;pho_not=deg2rad(15);
TP = 1/f;t0 = data(:, 3);t1 = t0 - t0(1);
t  = t1 ./ TP;
ycl     = data(:, 5) ; 
Ucl_mag = data(:, 6) ;Ucl_u   = data(:, 7);Ucl_v   = data(:, 8) ;
theta_var = data(:, 9);
% % % % % ------------------------------------------------------------
% % % % Normalisation 
% % % % ------------------------------------------------------------
NDU_v = (Ucl_v - mean(Ucl_v))/(h*f);
ND_ycl = (ycl - mean(ycl))/h;
NORM_theta = (theta_var - mean(theta_var))/pho_not;
theta_eq = mean(theta_var);
% % % % % ------------------------------------------------------------
% % % Plot: Phase-normalized signals
% % % ------------------------------------------------------------
figure(1);asetWH(gcf, 8, 8);legend off; hold on
fin = 199;
% % % % --- Left Y-axis ---
subplot(3,1,1);hold on
plot(t(1:fin), NORM_theta(1:fin), 'k:', 'LineWidth', lw);
plot(t(1:6:fin), NORM_theta(1:6:fin), 'bo', 'MarkerEdgeColor', 'k', ...
     'MarkerFaceColor', 'w', 'MarkerSize', ms);
ylabel('$\phi$', 'Interpreter', 'latex');
axis tight
xlabel('t/T');
asetWH(gcf, 6, 3);legend('','$f=5$','','f=100');
grid on; hold off
subplot(3,1,2);hold on
plot(t(1:fin), ND_ycl(1:fin), 'k:', 'LineWidth', lw);
plot(t(1:6:fin), ND_ycl(1:6:fin), 'bo', 'MarkerEdgeColor', 'b', ...
     'MarkerFaceColor', 'w', 'MarkerSize', ms);
ylabel('$\hat{y_{\mathrm{cl}}}$', 'Interpreter', 'latex');
axis tight
xlabel('t/T');
asetWH(gcf, 6, 6);legend off
grid on; hold off
subplot(3,1,3);hold on
plot(t(1:fin), NDU_v(1:fin), 'k:', 'LineWidth', lw);
plot(t(1:6:fin), NDU_v(1:6:fin), 'ro', 'MarkerEdgeColor', 'r', ...
     'MarkerFaceColor', 'w', 'MarkerSize', ms);
ylabel('$\hat{U_{\mathrm{cl}}}$', 'Interpreter', 'latex');
axis tight
xlabel('t/T');
asetWH(gcf, 9, 9);legend off
grid on; hold off


%% Plotting lissijous curves for position angle and velocity
clc;clear;close all
%%+++++++++++++++++++++++dt-12
fpath= "./data/D12_CL_info1.txt";
f = 300;h=9.13E-04;phi_not=deg2rad(10);
plot_hyst_loop(fpath,f,h,phi_not,'red')

%%+++++++++++++++++++++++dt-14
fpath= "./data/D14_CL_info1.txt";
f = 60;h=9.14E-05;phi_not=deg2rad(1);
plot_hyst_loop(fpath,f,h,phi_not,'blue')

function plot_hyst_loop(fpath,freq,h,phi_not,i) 
data = readmatrix(fpath,'NumHeaderLines',1); 
% data = readmatrix("../code_runs_Stokes/dt-2_d3.5f100/map_files/CL_info.txt"...
%     ,'NumHeaderLines',1); 
%% ------------------------------------------------------------
%  Data Processing
% ------------------------------------------------------------
f = freq;
TP = 1/f;
t0 = data(:, 3);        % 
t1 = t0 - t0(1);
t  = t1 ./ TP;
ycl     = data(:, 5)*1e3 ; % units in mm
Ucl_mag = data(:, 6);
Ucl_u   = data(:, 7) ;
Ucl_v   = data(:, 8) ;
theta_var = data(:, 9);

%% ------------------------------------------------------------
% Normalisation (same logic as Python)
% ------------------------------------------------------------
% NDU_v = (Ucl_v - mean(Ucl_v))/(h*f);
% NDU_ycl = (ycl - mean(ycl))/h ;
% NORM_theta = (theta_var - mean(theta_var))/phi_not;

NDU_v = (Ucl_v - mean(Ucl_v));
NDU_ycl = (ycl - mean(ycl)) ;
NORM_theta = (theta_var - mean(theta_var));

theta_eq = mean(theta_var);

%% %Normalized quantities
NORM_ycl = (ycl - mean(ycl));
U_mag = sign(Ucl_u) .* sqrt(Ucl_u.^2 + Ucl_v.^2);
NDU_mag = (U_mag - mean(U_mag)) ;
ND_theta = theta_var / abs(max(theta_var));
NORM_theta = (theta_var - mean(theta_var)) ;
NORM_v = (Ucl_v - mean(Ucl_v));


%% ---------------- Data sampling ----------------
endin = 199;
stp   = 5;

displacement = NDU_ycl(1:stp:endin);
velocity     = NORM_v(1:stp:endin);
angle        = NORM_theta(1:stp:endin);

ms = 8;   % marker size for overlays

%% ---------------- Figure handling ----------------
figure(1);asetWH(gcf,6,6);
set(gcf,'Color','w')
hold on

%% ================= (a) 3D Lissajous loop =================

% Backbone curve
plot3( velocity, displacement, angle, ...
      '-', 'Color', i, 'LineWidth',1.3)

% Markers
plot3( velocity, displacement, angle, ...
      'o', 'Color', i, 'MarkerSize', ms, ...
      'MarkerFaceColor', i,'MarkerEdgeColor','k')

xlabel('$U_{cl}$','Interpreter','latex')
ylabel('$\hat{y_{\mathrm{cl}}}$(mm)','Interpreter','latex')
zlabel('$\phi$','Interpreter','latex')

view(125,25)
grid('on')
box('on')
axis('tight')

% Apply 3D settings once
camproj('perspective')
camlight('headlight')
lighting('gouraud')
%% ================= (a) Velocity vs Angle =================
figure(2); hold on
ax4 = subplot(3,1,1);
hold(ax4,'on')

plot(ax4, displacement, angle, ...
     '-o', 'Color', i, 'LineWidth',1.3,'MarkerSize', ms, ...
     'MarkerFaceColor', i);


xlabel(ax4,'$\hat{y_{\mathrm{cl}}}$(mm)','Interpreter','latex')
ylabel(ax4,'$\phi$','Interpreter','latex')

axis(ax4,'tight')
grid(ax4,'on')
box(ax4,'on')


%% ================= (b) Velocity vs Displacement =================
ax2 = subplot(3,1,2);
hold(ax2,'on')

plot(ax2, velocity, displacement, ...
     '-o', 'Color', i, 'LineWidth',1.3,'MarkerSize', ms, ...
     'MarkerFaceColor', i)


xlabel(ax2,'$U_{cl}$','Interpreter','latex')
ylabel(ax2,'$\hat{y_{\mathrm{cl}}}$(mm)','Interpreter','latex')

axis(ax2,'tight')
grid(ax2,'on')
box(ax2,'on')


%% ================= (c) Velocity vs Angle =================
ax3 = subplot(3,1,3);
hold(ax3,'on')

plot(ax3, velocity, angle, ...
     '-o', 'Color', i, 'LineWidth',1.3, 'MarkerSize', ms, ...
     'MarkerFaceColor', i)


xlabel(ax3,'$U_{cl}$','Interpreter','latex')
ylabel(ax3,'$\phi$','Interpreter','latex')

axis(ax3,'tight')
grid(ax3,'on')
box(ax3,'on')



%######################################################################################
% --- Projections ---
% plot3(velocity, displacement, 0.8*min(angle)*ones(size(angle)), 'k--');  % XY plane
% plot3(velocity, 1*min(displacement)*ones(size(displacement)), angle, 'k--');  % XZ plane
% plot3(1.2*min(velocity)*ones(size(velocity)), displacement, angle, 'k--');  % YZ plane
hold off;
end
figure(1);asetWH(gcf,6,6);legend('','Split','','Rolling');
figure(2);asetWH(gcf,9,9);legend off
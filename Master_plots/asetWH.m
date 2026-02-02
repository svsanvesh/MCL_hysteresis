function  asetWH(h, width,height )
pos = get(h, 'Position');
set(h, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size
% Here we preserve the size of the image when we save it.
set(h,'InvertHardcopy','on');
set(h,'PaperUnits', 'inches');
papersize = get(h, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(h,'PaperPosition', myfiguresize);
set(findall(gcf,'-property','FontSize'),'FontSize',18,'FontWeight','bold') % adjust fontsize to your document
set(findall(gcf,'-property','Box'),'Box','on') % optional
set(findall(gcf,'-property','Interpreter'),'Interpreter','latex') 
set(findall(gcf,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(gca, 'TickLength', [0.02, 0.01]); % Major tick: 0.02, Minor tick: 0.01
set(gca, 'FontSize', 18,'FontWeight','bold'); % Adjust the value as needed
set(gca, 'LineWidth', 1.3);
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on'); % Turn on minor ticks for both axes
% set(gca, 'MinorGridLineStyle', '--'); % Optional: customize the style of minor grid lines
hLegend = legend; % Get the legend handle
set(hLegend, 'LineWidth', 0.8); % Adjust the line width of the legend box
%% %%%%%% Construct the new filename
% print('The plot is to be saved as a pdf and .fig file')
% fname = input('enter pdf name:');
% figfname = input('enter fig name:');
% savefig(figfname);
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print(gcf, '-dpdf','-r600', fname);
end

% print(gcf,'-dpng','-r300','test.png')


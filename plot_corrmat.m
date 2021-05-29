function [plot_s] = plot_corrmat(Cor,Name)
% PLOT_CORRMAT  % Calculate and/or Plot Correlation Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Produce the input lower triangular matrix data:
C=tril(Cor);
%% Compute center of each circle
%% This assumes the x and y values were not entered in imagesc()
%% Set color of each rectangle
%% Set color scale
%#cmap = parula(256);%jet
%#cmap(1,:) = 1;
%#cmap(end,:) = 1;
%#colormap(cmap);
cmap = colormap (summer (32));
imagesc(C); % plot the matrix
caxis([0, 1]);
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
% warning ("off")
set(gcf, 'Toolbar', 'none', 'Menu', 'none');
myLabel = cellstr(['--';'S1';'--';'S2';'--';'S3';'--']);
set(gca, 'XTickLabel', myLabel);
set(gca, 'YTickLabel', myLabel);
% xticklabels({'S1','S2','S3'});
% yticklabels({'S1','S2','S3'});
% yticks(['S1','S2','S3']);
% xticks(['S1','S2','S3']);
colorbar;
title(Name);
% cmap = CorrColormap; % Uncomment for CorrColormap
% Cscaled = (C - clrLim(1))/range(clrLim); % always [0:1]
% colIdx = discretize(Cscaled,linspace(0,1,size(cmap,1)));
% Set size of each circle
% Scale the size between [0 1]
% Cscaled = (abs(C) - 0)/1;
% Create figure
% ax = axes(fh);
% ax = axes(fh,'Position',[0.5 0.5 3.5 3.5],'Box','on');
% hold(ax,'on')
% colormap(ax,'jet');
% colormap(CorrColormap) %Uncomment for CorrColormap
% diamLim = [0.1, 1];
% diamSize = Cscaled * range(diamLim) + diamLim(1);
% tickvalues = 1:length(C);
% x = zeros(size(tickvalues));
% text(x, tickvalues, myLabel, 'HorizontalAlignment', 'right');
% x(:) = length(C)+1;
% text(tickvalues, x, myLabel, 'HorizontalAlignment', 'right','Rotation',90);
% % Create circles
% theta = linspace(0,2*pi,25); % the smaller, the less memory req'd.
% h = arrayfun(@(i)fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
%     diamSize(i)/2 * sin(theta) + yAll(i), cmap(colIdx(i),:),'LineStyle','none'),1:numel(xAll));
% axis(ax,'equal')
% axis(ax,'tight')
% set(ax,'YDir','Reverse')
% colorbar()
% caxis(clrLim);
% imagesc(C)
% axis off
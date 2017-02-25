function text_handles = featureticks(handle) 
% Noam Cohen,
% Based on polarticks function by Adam Danz
fSize = 22;
if nargin < 3 || isempty(handle)
    handle = gca; end
    
% remove junck
h = findall(handle,'type','line');          %get handles for all lines in polar plot
delete(h((1:12)));

% delete text
t = findall(handle,'type','text');          %get handles for all text in polar plot
delete (t)

v = [get(handle, 'XLim') get(handle, 'YLim')];
rmax = v(2);
rt = 1 * rmax./11.6906729634002;

% Create Gray text
text('Parent',handle,'FontWeight','bold','FontSize',fSize,'String',' Off',...
    'Position',rt.*[11.6906729634002 3.87012987012987 0],...
    'Visible','on');

% Create Chirp text
text('Parent',handle,'HorizontalAlignment','right','FontWeight','bold','FontSize',fSize,'String','Chirp',...
    'Position',rt.*[-6.91243604879969 -10.0958284140102 0],...
    'Visible','on');

% Create On text
text('Parent',handle,'VerticalAlignment','bottom','HorizontalAlignment','center','FontWeight','bold',...
    'FontSize',fSize,'String','Amp','Position',rt.*[-12.6906729634002 3.87012987012987 0],'Visible','on');

% Create Off text
text('Parent',handle,'VerticalAlignment','top','HorizontalAlignment','center','FontWeight','bold','FontSize',fSize,...
    'String','On','Position',rt.*[0.332152695789058 13.5367965367965 0],'Visible','on');

% Create Amp text
text('Parent',handle,'FontWeight','bold','FontSize',fSize,'String','Gray',...
    'Position',rt.*[6.91243604879969 -10.0958284140102 0],...
    'Visible','on');

% Create title
title('Chirp Feature Visual','FontWeight','bold');

% Create colorbar
colorbar('peer',handle,'Position',...
    [0.862415730337079 0.102877350912666 0.0467977528089888 0.334249231365815],...
    'Ticks',[0 0.0833333333333333 0.166666666666667 0.25]);

end

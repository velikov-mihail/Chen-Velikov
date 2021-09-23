function setplot(titletext,xtext,ytext,scale,font,interpreter)

% setplot() or 
% setplot(scale)

if nargin <= 0
    titletext = '';
end
if nargin <= 1
    xtext = '';
end
if nargin <= 2
    ytext = '';
end
if nargin <= 3
    scale = 1;
end
if nargin <= 4    
    font = 'Helvetica';
end
if nargin <= 5
    interpreter = 'none';
end



%%
if nargin > 0
    title(titletext)
    xlabel(xtext)
    ylabel(ytext)
end

%%



set(gcf, 'Color', 'w');

% Tick labels: must be done first
set(gca,'fontsize', scale*10);


% Title 
tithand = get(gca,'title');
set(tithand, ...
    'FontName'   , font,...
    'FontSize'   , scale*13,...
    'FontWeight' , 'normal', ...
    'interpreter', interpreter);


% Axis Labels
xlhand = get(gca,'xlabel');
ylhand = get(gca,'ylabel');

set(xlhand, ...
    'FontName'   , font,...
    'FontSize'   , scale*12 ,...
    'interpreter', interpreter);

set(ylhand, ...
    'FontName'   , font,...
    'FontSize'   , scale*12 ,...
    'interpreter', interpreter);



% legend
hlegend = findobj(gcf,'Type','axes','Tag','legend');
set(hlegend, ...
    'FontName'   , font,...
    'FontSize'   , scale*12 ,...
    'interpreter', interpreter);





%% lines
h = get(gca, 'children');
for i = 1:length(h)
    temp = get(h(i));
    if isfield(temp,'XData') & ~isfield(temp,'BarWidth')
        set(h(i), 'linewidth'  , scale*1, ...
                'markersize', scale*6);
    end
end        

%%

% Axis lines
set(gca,'linewidth', scale*1.5);






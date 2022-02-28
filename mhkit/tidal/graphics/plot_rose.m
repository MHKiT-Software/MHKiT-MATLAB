function figure=plot_rose(data, width_dir, width_vel, ...
                                        options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Creates a polar histogram. Direction angles from binned histogram must 
%       be specified such that 0  degrees is north.
%     
% Parameters
% ------------  
%    data: structure
%
%      data.time: vector
%       days from January 0, 0000 in the proleptic ISO calendar
%
%      data.d: vector
%       time-series of directions [degrees]
%
%      data.s: vector
%       time-series of speeds [cm/s]
%
%    width_dir: float 
%        Width of directional bins for histogram in degrees
%
%    width_vel: float 
%        Width of velocity bins for histogram in m/s
%
%    flood_ebb: 2 element vector (optional)
%        Direction in degrees added to theta ticks
%        to call: plot_rose(Q, width_dir, width_vel,"flood_ebb",flood_ebb)
%
%    title: string (optional)
%        title for the plot 
%        to call: plot_rose(Q, width_dir, width_vel,"title",title)
%
%    savepath: string (optional)
%        path and filename to save figure.
%        to call: plot_rose(Q, width_dir, width_vel,"savepath",savepath)
%
% 
% Returns
% ---------
%   figure handle to water current rose plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    data
    width_dir
    width_vel
    options.title = "";
    options.savepath = "";
    options.flood_ebb = [];
end

% checking input variables
% if nargin > 5
%     ME = MException('MATLAB:plot_rose','Incorrect number of input arguments, too many agruments, requires 5 at most, %d arguments passed',nargin);
%     throw(ME);
% end;

%check to see if the first input argumeent is a structure
if any(~isstruct(data))
    ME = MException('MATLAB:plot_rose','data must be a structure');
    throw(ME);
end

%check to see if the second input argumeent is a number
if any([~isnumeric(width_dir), length(width_dir) ~= 1])
    ME = MException('MATLAB:plot_rose','width_dir must be a number');
    throw(ME);
end

%check to see if the second input argumeent is a number
if any([~isnumeric(width_vel), length(width_vel) ~= 1])
    ME = MException('MATLAB:plot_rose','width_vel must be a number');
    throw(ME);
end

% Parsing out the variable arguments
% if nargin >= 4
%     for indx = 4:1:nargin
%         if (ischar(varargin{indx-3}) | isstring(varargin{indx-3}))
%             titleText = varargin{indx-3};
%         elseif isnumeric(varargin{indx-3})
%             if length(varargin{indx-3}) ==2
%                 flood_ebb = varargin{indx-3};
%             end;
%         end;
%     end;
% end;
    
% Plot Variables
radialLabelDir = 30;
propRingRad    = 5;
CMaPType = 'spring';

% Calculating the 2D histogram
DirEdges = 0:width_dir:360;
VelEdges = 0:width_vel:max(data.s);

% converting degrees to radians
theta = data.d/180*pi;

% setting the colors
cmap = colormap(CMaPType);
colors = interp1(linspace(0,max(data.s),length(cmap)),cmap,VelEdges);

%tx = Q.s'<Yedges(end)

figure = polarhistogram(theta(data.s<VelEdges(end)),DirEdges/180*pi,'FaceColor',colors(end,:), ...
        'displayname',[num2str(VelEdges(end-1)) ' - ' num2str(VelEdges(end)) 'm/s']);
hold on
for idx = length(VelEdges)-1:-1:2
    polarhistogram(theta(data.s<VelEdges(idx)),DirEdges/180*pi,'FaceColor',colors(idx,:), ...
            'displayname',[num2str(VelEdges(idx-1)) ' - ' num2str(VelEdges(idx)) '\itm/s']);
end

% Adding legend to plot
legend('Show');

% adding title
title(options.title)

[N] = histcounts2(data.d,data.s,DirEdges,VelEdges,'Normalization','probability');
%forming the vector with the total probability for each direction bin
dir_prob = sum(N,2)*100;

% Setup the axes of the polar plot
axinfo = gca;
axinfo.ThetaZeroLocation = 'top';
axinfo.ThetaDir = 'clockwise';
axinfo.ThetaTick = 0:45:360;
axinfo.ThetaTickLabel = {'N','NE','E','SE','S','SW','W','NW'};

axinfo.RAxisLocation = radialLabelDir;
RadRingLoc = propRingRad:propRingRad:propRingRad*ceil(max(dir_prob)/propRingRad);
axinfo.RTick = RadRingLoc/100*length(data.s);
axinfo.RTickLabel = num2cell(RadRingLoc);
axinfo.FontSize = 20;

% adding in the flood and ebb directions if they are included

if length(options.flood_ebb) == 2
    % adding in the radial lines for the ebb and flood 
    polarplot([1,1]*options.flood_ebb(1)/180*pi,[0,max(RadRingLoc/100)*length(data.s)],'linewidth',2,'color','r','displayname','Flood Direction');
    text(options.flood_ebb(1)/180*pi,max(RadRingLoc/100)*length(data.s),'Flood','fontsize',15,'HorizontalAlignment','Center','BackGroundColor','w');
    polarplot([1,1]*options.flood_ebb(2)/180*pi,[0,max(RadRingLoc/100)*length(data.s)],'linewidth',2,'color','b','displayname','Ebb Direction');    
    text( options.flood_ebb(2)/180*pi,max(RadRingLoc/100)*length(data.s),'Ebb','fontsize',15,'HorizontalAlignment','Center','BackGroundColor','w');
end

len = strlength(options.savepath);
if len > 1
    saveas(figure, options.savepath);
end

hold off   


% if nargin == 4
%     title(varargin{1})
% end
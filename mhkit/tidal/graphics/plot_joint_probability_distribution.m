function fig=plot_joint_probability_distribution(Q, width_dir, width_vel, ...
                                        options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Creates a polar histogram. Direction angles from binned histogram must 
%   be specified such that 0 is north.
%     
% Parameters
% ------------   
%    data: structure
%
%     data.time: vector
%           days from January 0, 0000 in the proleptic ISO calendar
%
%     data.d: vector
%           time-series of directions [degrees]
%
%     data.s: vector
%           time-series of speeds [cm/s]
%
%    width_dir: float 
%       Width of directional bins for histogram in degrees
%
%    width_vel: float 
%       Width of velocity bins for histogram in m/s
%
%    flood_ebb: 2 element vector (optional)
%       Direction in degrees added to theta ticks 
%       to call: plot_joint_probability_distribution(Q, width_dir, width_vel,"flood_ebb",flood_ebb)
%
%    title: string (optional)
%       title for the plot 
%       to call: plot_joint_probability_distribution(Q, width_dir, width_vel,"title",title)
%
%    savepath: string (optional)
%       path and filename to save figure.
%       to call: plot_joint_probability_distribution(Q, width_dir, width_vel,"savepath",savepath)
% 
% Returns
% ---------
%   figure handle to joint probability distribution plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    Q
    width_dir
    width_vel
    options.flood_ebb = [];
    options.title = "";
    options.savepath = "";

end

% checking input variables
% if nargin > 5
%     ME = MException('MATLAB:plot_joint_probability_distribution','Incorrect number of input arguments, too many agruments, requires 5 at most, %d arguments passed',nargin);
%     throw(ME);
% end;

%check to see if the first input argumeent is a structure
if any(~isstruct(Q))
    ME = MException('MATLAB:plot_joint_probability_distribution','Q must be a structure');
    throw(ME);
end;

%check to see if the second input argumeent is a number
if any([~isnumeric(width_dir), length(width_dir) ~= 1])
    ME = MException('MATLAB:plot_joint_probability_distribution','width_dir must be a number');
    throw(ME);
end;

%check to see if the second input argumeent is a number
if any([~isnumeric(width_vel), length(width_vel) ~= 1])
    ME = MException('MATLAB:plot_joint_probability_distribution','width_vel must be a number');
    throw(ME);
end;

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

% matlab does not have the ability to overlay a pcolor on polar axes, so a
% polar plot must manually be created.

% Plot Variables
max_rings = 8; % maximum number of anular lines
Vel_label_ang = 80/180*pi; %direction of the velocity lables
pcolorThreshold = 0.05; %threshold for the smallest number to include in pcolor
CMaPType = 'jet';

% converting degrees to radians
theta = Q.d/180*pi;

% Calculating the 2D probability
DirEdges = 0:width_dir:360-width_dir;
VelEdges = 0:width_vel:max(Q.s);
[H] = histcounts2(Q.d,Q.s,DirEdges,VelEdges,'Normalization','probability');
H = H*100;

% creating vel & dir bins index to point to the middle of bin 
dir_bins = DirEdges(2:end)-width_dir/2;
dir_bins(1) = 0;
vel_bins = VelEdges(2:end)-width_vel/2;
vel_bins(1) = 0;

% converting data to cartesian
Xs = Q.s.*sin(theta);
Ys = Q.s.*cos(theta);

% plotting the scatterpoints
fig = plot(Xs,Ys,'.k');
hold on
set(gca,'xtick',[],'ytick',[]);

% setting up the axes
% circles
MaxR = round(max(Q.s),1);
rectangle('Position',[-MaxR -MaxR MaxR*2 MaxR*2],'Curvature',[1,1]);
idRng = ceil(length(VelEdges)/max_rings);
for idx = 2:idRng:length(VelEdges)-1
    rectangle('Position',[-VelEdges(idx) -VelEdges(idx) VelEdges(idx)*2 VelEdges(idx)*2],'Curvature',[1,1]);
    text(VelEdges(idx)*sin(Vel_label_ang),VelEdges(idx)*cos(Vel_label_ang),[num2str(VelEdges(idx)) '\itm/s'],'fontsize',12,'HorizontalAlignment','Center','BackGroundColor','w');
end;

% radial lines and axis labels
Xrad = MaxR.*sin((0:45:359)/180*pi);
Yrad = MaxR.*cos((0:45:359)/180*pi);
AxisLabel = {'N','NE','E','SE','S','SW','W','NW'};
for idx = 1:8
    line([0 Xrad(idx)],[0 Yrad(idx)]);
    text(Xrad(idx)*1.1,Yrad(idx)*1.05,AxisLabel(idx),'fontsize',15,'HorizontalAlignment','Center');
end

% adding title
title(options.title)

% adding in the flood and ebb directions if they are included
if length(options.flood_ebb) == 2 
    line([0,MaxR*sin(options.flood_ebb(1)/180*pi)],[0 MaxR*cos(options.flood_ebb(1)/180*pi)]);
    text(MaxR*sin(options.flood_ebb(1)/180*pi)*1.1,MaxR*cos(options.flood_ebb(1)/180*pi)*1.1,'Flood','fontsize',15,'HorizontalAlignment','Center','BackGroundColor','w');
    line([0,MaxR*sin(options.flood_ebb(2)/180*pi)],[0 MaxR*cos(options.flood_ebb(2)/180*pi)]);
    text(MaxR*sin(options.flood_ebb(2)/180*pi)*1.1,MaxR*cos(options.flood_ebb(2)/180*pi)*1.1,'Ebb','fontsize',15,'HorizontalAlignment','Center','BackGroundColor','w');
end


% creating direction and velocity matrix as indexes into the probability
% matrix
dir_matrix = ones(size(H)).*(dir_bins'/180*pi);
y_pos = cos(dir_matrix).*vel_bins;
y_pos = [y_pos; y_pos(1,:)]; 
x_pos = sin(dir_matrix).*vel_bins;
x_pos = [x_pos; x_pos(1,:)]; 
H = [H;H(1,:)];

% creating the color plot for the PDF
Hplt = H;
Hplt(H<pcolorThreshold) = nan;
pfig = pcolor(x_pos,y_pos,Hplt);
set(pfig, 'EdgeColor', 'none');
colormap(CMaPType)
shading interp;
cbh = colorbar;
ylabel(cbh,'Joint Probability [%]','fontsize',20);

len = strlength(options.savepath);
if len > 1
    saveas(figure, options.savepath);
end 

hold off


% if nargin == 4
%     title(varargin{1})
% end
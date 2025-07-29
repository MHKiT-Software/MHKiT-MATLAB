function plot_mooring_animate(data, dimension, step, varargin)
%Graphics function that creates a 2D or 3D animation of the node 
% positions of a mooring line over time.
%   Detailed explanation goes here

arguments (Input)
    data
    dimension {mustBeTextScalar} = '2D'
    step {mustBeNumeric} = 0.05
end
arguments (Repeating)
    varargin
end

%channel names
chans = data.Properties.VariableNames;

% get node names
idx = strfind(chans, 'x');
idx = cell2mat(cellfun(@(x)any(~isempty(x)),idx,'UniformOutput',false));
nodes_x = {chans{idx}};
idy = strfind(chans, 'y');
idy = cell2mat(cellfun(@(x)any(~isempty(x)),idy,'UniformOutput',false));
nodes_y = {chans{idy}};
idz = strfind(chans, 'z');
idz = cell2mat(cellfun(@(x)any(~isempty(x)),idz,'UniformOutput',false));
nodes_z = {chans{idz}};

figure;
% allow for axis labels
for k = 1:2:length(varargin)
    labelType = varargin{k};
    labelText = varargin{k+1};
    
    switch lower(labelType)
        case 'title'
            title(labelText);
        case 'xlabel'
            xlabel(labelText);
        case 'ylabel'
            ylabel(labelText);
        case 'zlabel'
            zlabel(labelText);
        otherwise
            warning('Unknown label type: %s', labelType);
    end
end

% 2D animate
if strcmp(dimension,'2D')
    h = animatedline('Marker','o');
    grid on
    clearpoints(h);
    for k=1:size(data,2)
        clearpoints(h);
        x = table2array(data(k, nodes_x));
        z = table2array(data(k, nodes_z));
        addpoints(h,x,z)
        drawnow
        pause(step)
    end
% 3D animate
elseif strcmp(dimension,'3D')
    view(3);
    h = animatedline('Marker','o');
    grid on
    clearpoints(h);
    for k=1:size(data,2)
        clearpoints(h);
        x = table2array(data(k, nodes_x));
        y = table2array(data(k, nodes_y));
        z = table2array(data(k, nodes_z));
        addpoints(h,x,y,z)
        drawnow
        pause(step)
    end
else
    error('Please make sure to specify the correct dimension!')
end



end


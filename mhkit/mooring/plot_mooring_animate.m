function plot_mooring_animate(data, dimension)
%Graphics function that creates a 2D or 3D animation of the node 
% positions of a mooring line over time.
%   Detailed explanation goes here

arguments
    data
    dimension {mustBeTextScalar} = '2d'
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

% 2D animate
if strcmp(dimension,'2d')
    h = animatedline;
    for k=1:size(data,2)
        clearpoints(h)
        x = table2array(data(k, nodes_x));
        z = table2array(data(k, nodes_z));
        addpoints(h,x,z)
        drawnow
    end
end

end


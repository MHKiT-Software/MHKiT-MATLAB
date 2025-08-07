function line_lay_length = lay_length(data, depth, tolerance)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the lay length of a mooring line over time.
%
% Lay length is the measure of how much of the mooring line is in contact 
% with the seabed. This function requires MoorDyn line output data containing
% individual node positions NodeNpx, NodeNpy, NodeNpz.
%
%   Parameters
%    ----------
%    data: table
%        MoorDyn line output data table containing Time column and node 
%        position data (NodeNpx, NodeNpy, NodeNpz) where N is the node number.
%        Data should be obtained using mhkit.mooring.read_moordyn() with a 
%        MoorDyn line output file (*.Line*.out).
%    depth: double
%        Depth of seabed (m). Should be negative for depths below sea level.
%    tolerance: double, optional
%        Tolerance to detect first lift point from seabed, by default 0.25
%        meters. Nodes with z-position > depth + abs(tolerance) are 
%        considered lifted from seabed.
%
%    Returns
%    -------
%    line_lay_length: double array
%        Array containing the lay length at each time step (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


arguments
    data table
    depth {mustBeNumeric}
    tolerance {mustBeNumeric} = 0.25
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

% check if data contains the necessary nodes
if isempty(nodes_x) || isempty(nodes_y) || isempty(nodes_z)
    error("The data must contain x, y, and z node data!")
end
if length(nodes_z) < 3
    error("This function requires at least 3 nodes to calculate lay length")
end

% find name of first z point where tolerance is exceeded
zdata = table2array(data(:,nodes_z)); % Convert table to array
laypoint_mask = zdata > (depth + abs(tolerance));
laypoint = nodes_z(any(laypoint_mask));

% get previous z-point
lay_indx = find(cellfun(@(x) strcmp(x, laypoint), nodes_z)) - 1;
lay_z = nodes_z(lay_indx);

% get corresponding x-point and y-point node names
lay_x = strcat(lay_z{:}(1:end-1), "x");
lay_y = strcat(lay_z{:}(1:end-1), "y");
lay_0x = nodes_x{1};
lay_0y = nodes_y{1};

% find distance between initial point and lay point
laylength_x = data.(lay_x) - data.(lay_0x);
laylength_y = data.(lay_y) - data.(lay_0y);
line_lay_length = (laylength_x.^2 + laylength_y.^2) .^ 0.5;

end



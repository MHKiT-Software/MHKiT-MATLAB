function line_lay_length = lay_length(data, depth, tolerance)
%Calculate the laylength of a mooring line over time.
%   Parameters
%    ----------
%    data: Array
%        Array containing x,y,z nodes (ie Node1px, Node1py, Node1pz)
%    depth: double
%        Depth of seabed (m)
%    tolerance: double, optional
%        Tolerance to detect first lift point from seabed, by default 0.25
%        meters
%
%    Returns
%    -------
%    line_lay_length: xr.Dataset
%        Array containing the laylength at each time step

arguments
    data {mustBeNumeric}
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
if isempty(nodes_x) || isempty(nodes_y) || isempty(nodes_zz)
    error("The data must contain x, y, and z node data!")
end
if length(nodes_z) < 3
    error("This function requires at least 3 nodes to calculate lay length")
end

% find name of first z point where tolerance is exceeded
laypoint = data(:,nodes_z) > depth + abs(tolerance);
laypoint = nodes_z(any(table2array(laypoint)));

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



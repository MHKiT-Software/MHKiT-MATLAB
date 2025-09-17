function fmin = minimum_frequency(water_depth, c, c_seabed)
% 
% Estimate the shallow water cutoff frequency based on the speed of
% sound in the water column and the speed of sound in the seabed
% material (generally ranges from 1450 - 1800 m/s)
% 
% Parameters
% ----------
% water_depth: double
%     Depth of the water column in meters.
% c: double, optional
%     Speed of sound in the water column in meters per second. Default is 1500 m/s.
% c_seabed: double, optional
%     Speed of sound in the seabed material in meters per second. Default is 1700 m/s.
% 
% Returns
% -------
% f_min: double
%     The minimum cutoff frequency in Hz.
% 
% Reference
% ---------
% Jennings 2011 - Computational Ocean Acoustics, 2nd ed.


arguments (Input)
    water_depth {mustBeNumeric}
    c {mustBeNumeric} = 1500
    c_seabed {mustBeNumeric} = 1700
end

arguments (Output)
    fmin {mustBeNumeric}
end

water_depth = double(water_depth);

if any(water_depth <= 0)
    error('All elements of water_depth must be positive numbers.');
end
if c <= 0
    error('c must be a positive number.');
end
if c_seabed <= 0
    error('c_seabed must be a positive number.');
end
if c_seabed <= c
    error('c_seabed must be greater than c.');
end

fmin = c ./ (4 .* water_depth .* sqrt(1 - (c ./ c_seabed).^2));

end


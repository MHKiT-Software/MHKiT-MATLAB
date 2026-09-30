function Fr = froude_number(v, h, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the Froude number of the river, channel or duct flow
%
% Used to check the subcritical flow assumption (Fr < 1).
%
% Parameters
% ------------
% v : double [m/s]
%   Average velocity
% h : double [m]
%   Mean hydraulic depth
% g : double [m/s^2] (optional)
%   Name-value argument. Gravitational acceleration, default 9.80665
%
% Returns
% ---------
% Fr : double [-]
%   Froude number of the river
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    v (1,1) {mustBeNumeric}
    h (1,1) {mustBeNumeric}
    options.g (1,1) {mustBeNumeric} = 9.80665
end

arguments (Output)
    Fr (1,1) {mustBeNumeric}
end

Fr = v / sqrt(options.g * h);

end

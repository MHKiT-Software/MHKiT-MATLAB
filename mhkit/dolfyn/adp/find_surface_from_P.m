function out = find_surface_from_P(ds, options)
% Calculates the distance to the water surface. Temperature and salinity
% are used to calculate seawater density, which is in turn used with the
% pressure data to calculate depth.
% 
% Parameters
% ----------
% ds: Dataset
%   The full adcp dataset
% salinity: numeric
%   Water salinity in psu
% 
% Returns
% -------
% out: Dataset
%   adds the variables "water_density" and "depth" to the input dataset.
% 
% Notes
% -----
% Requires that the instrument's pressure sensor was calibrated/zeroed
% before deployment to remove atmospheric pressure.
% 
% Calculates seawater density at normal atmospheric pressure according
% to the UNESCO 1981 equation of state. Does not include hydrostatic 
% pressure.

    arguments
        ds         
        options.salinity = 35;
    end

    % Density calcation
    T = ds.temp.data;
    S = options.salinity;
    % standard mean ocean water
    rho_smow = 999.842594 + 6.793953e-2*T - 9.095290e-3*T.^2 + ...
        1.001685e-4*T.^3 - 1.120083e-6*T.^4 + 6.536332e-9*T.^5;
    % at water surface
    B1 = 8.2449e-1 - 4.0899e-3*T + 7.6438e-5*T.^2 - 8.2467e-7*T.^3 + ...
        5.3875e-9*T.^4;
    C1 = -5.7246e-3 + 1.0227e-4*T - 1.6546e-6*T.^2;
    d0 = 4.8314e-4;
    rho_atm0 = rho_smow + B1*S + C1*S^1.5 + d0*S^2;

    % Depth = pressure (conversion from dbar to MPa) / water weight
    d = (ds.pressure.data*10000)./(9.81*rho_atm0);

    if isfield(ds.attrs, 'h_deploy')
        d = d + ds.attrs.h_deploy;
        description = "Water depth to seafloor";
    else
        description = "Water depth to ADCP";
    end

    ds.water_density.data = rho_atm0;
    ds.water_density.dims = { 'time' };
    ds.water_density.coords.time = ds.coords.time;
    ds.water_density.units = "kg/m^3";
    ds.water_density.description = ['Water density according to UNESCO '...
        '1981 equation of state'];

    ds.depth.data = d;
    ds.depth.dims = { 'time' };
    ds.depth.coords.time = ds.coords.time;
    ds.depth.units = "m";
    ds.depth.description = description;

    out = ds;

end


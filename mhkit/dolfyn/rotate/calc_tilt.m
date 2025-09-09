function tilt = calc_tilt(pitch, roll, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate instrument tilt from pitch and roll angles.
%
% This function calculates the total tilt angle of an ADCP instrument from
% pitch and roll measurements. This function warns the user if the tilt is above 5°
% as tilts above this threshold are likely to have a negative affect on the accuracy 
% of flow and turbulence calculations.
%
% Parameters
% ------------
%   pitch: array
%       time series of pitch angle (forward/backward tilt)
%   roll: array  
%       time series of roll angle (side-to-side tilt)
%   units: string
%       Units of input angles: 'degrees', 'deg', 'radians' or 'rad'
%       Default: 'degrees'
%   output_units: string
%       Units for output tilt: 'degrees', 'deg', 'radians' or 'rad'
%       Default: same as input units
%
% Returns
% ---------
%   tilt: array
%       tilt angle in specified output units
%
% Algorithm
% ---------
% tilt_rad = atan( √( tan(roll_rad)² + tan(pitch_rad)² ) )
%
% Example
% -------
% % Calculate a time series of tilt from pitch and roll time series in degrees
% tilt = calc_tilt(pitch_data, roll_data, 'units', 'degrees');
%
% Notes
% -----
% - Large tilts (> 5°) can affect turbulence measurements
% - This function issues warnings for tilts > 5°
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    arguments
        pitch
        roll
        options.units = 'degrees'
        options.output_units = ''
    end
    
    % Validate inputs
    if ~isnumeric(pitch) || ~isnumeric(roll)
        error('mhkit:dolfyn:calc_tilt: Pitch and roll must be numeric arrays');
    end
    
    if ~isequal(size(pitch), size(roll))
        error('mhkit:dolfyn:calc_tilt: Pitch and roll arrays must have the same size');
    end
    
    % Validate units
    valid_units = {'degrees', 'radians', 'deg', 'rad'};
    if ~ismember(lower(options.units), valid_units)
        error('mhkit:dolfyn:calc_tilt: Units must be ''degrees'' or ''radians''');
    end
    
    % Set default output units
    if isempty(options.output_units)
        options.output_units = options.units;
    end
    
    if ~ismember(lower(options.output_units), valid_units)
        error('mhkit:dolfyn:calc_tilt: Output units must be ''degrees'' or ''radians''');
    end
    
    % Normalize unit names
    input_is_degrees = ismember(lower(options.units), {'degrees', 'deg'});
    output_is_degrees = ismember(lower(options.output_units), {'degrees', 'deg'});
    
    % Convert to radians if needed for calculation
    if input_is_degrees
        pitch_rad = deg2rad(pitch);
        roll_rad = deg2rad(roll);
    else
        pitch_rad = pitch;
        roll_rad = roll;
    end
    
    % Calculate tilt using trigonometric relationship
    % Source: mhkit_python:dolfyn.rotate.base.calc_tilt, author @jmcvey3
    tilt_rad = atan(sqrt(tan(roll_rad).^2 + tan(pitch_rad).^2));

    % Quality assessment and warnings
    % Convert tilt to degrees for quality assessment regardless of output units
    tilt_deg = rad2deg(tilt_rad);
    max_tilt_deg = max(tilt_deg(:));
    mean_tilt_deg = mean(tilt_deg(:));
    
    % Issue warnings for large tilts
    if max_tilt_deg > 5
        warning('mhkit:dolfyn: Maximum tilt %.1f° exceeds recommended 5° limit for accurate turbulence measurements', max_tilt_deg);
        fprintf('Tilt Analysis Summary:\n');
        fprintf('  Mean tilt: %.2f°\n', mean_tilt_deg);
        fprintf('  Max tilt:  %.2f°\n', max_tilt_deg);
        fprintf('  Std tilt:  %.2f°\n', std(tilt_deg(:)));
    end
    
    % Convert to desired output units
    if output_is_degrees
        tilt = rad2deg(tilt_rad);
    else
        tilt = tilt_rad;
    end
end

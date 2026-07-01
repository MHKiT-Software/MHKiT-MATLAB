function spsd_cal = apply_calibration(spsd,sensitivity_curve, fill_value, interp_method)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Apply custom calibration to spectral density values
%
% Parameters
% ------------
%   spsd: struct
%       Mean square sound pressure spectral density in V^2/Hz
%       spsd.data : Spectral density data [V^2/Hz]
%       spsd.freq : Frequency vector [Hz]
%       spsd.time : Time vector (if time series)
%   sensitivity_curve: matrix
%       Calibrated sensitivity curve in units of dB rel 1 V^2/uPa^2
%       First column should be frequency, second column should be calibration values
%   fill_value: float
%       Value with which to fill missing values from the calibration curve [dB rel 1 V^2/uPa^2]
%   interp_method: string
%       Interpolation method to use when interpolating the calibration
%       curve to the frequencies in 'spsd'. Default is "linear". 
%
% Returns
% ---------
%   spsd_cal: struct
%       Calibrated spectral density in Pa^2/Hz, indexed by time and frequency
%       spsd_cal.data : Calibrated spectral density data [Pa^2/Hz]
%       spsd_cal.freq : Frequency vector [Hz]
%       spsd_cal.time : Time vector (if time series)
%       spsd_cal.name : Data name
%       spsd_cal.units : Data units
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    spsd struct
    sensitivity_curve {mustBeNumeric, mustBeFinite}
    fill_value {mustBeNumeric, mustBeFinite}
    interp_method {mustBeTextScalar} = "linear"
end

arguments (Output)
    spsd_cal struct
end

% Validate spsd structure
validate_spsd_struct(spsd, 'apply_calibration', ...
    'required_fields', {{'freq'}});

spsd_cal = spsd;

% interpolate calibration
calibration = interp1(sensitivity_curve(:,1), ...
    sensitivity_curve(:,2), spsd.freq, interp_method);
% use first cal value to fill NaNs in lower freq
idx = find(spsd.freq < sensitivity_curve(1,1));
calibration(idx) = fillmissing(calibration(idx),'constant',sensitivity_curve(1,2));
% fill NaNs at high freq using specified fill_value
calibration = fillmissing(calibration,'constant',fill_value);

% subtract from sound pressure spectral density
sense_ratio = 10.^(calibration / 10); % V^2/uPa^2
spsd_cal.data = (spsd_cal.data ./ sense_ratio) / 1e12; % Pa^2/Hz
spsd_cal.name = "Calibrated Sound Pressure Spectral Density";
spsd_cal.units = "Pa^2/Hz";

end

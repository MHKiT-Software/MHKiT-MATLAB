function spsd_cal = apply_calibration(spsd,sensitivity_curve, fill_value)

% Applies custom calibration to spectral density values.
% 
% Parameters
% ----------
% spsd: struct
%     Mean square sound pressure spectral density in V^2/Hz.
% sensitivity_curve: matrix
%     Calibrated sensitivity curve in units of dB rel 1 V^2/uPa^2.
%     First column should be frequency, second column should be calibration values.
% fill_value: float or int
%     Value with which to fill missing values from the calibration curve,
%     in units of dB rel 1 V^2/uPa^2.
% 
% Returns
% -------
% spsd_calibrated: struct
%     Spectral density in Pa^2/Hz, indexed by time and frequency.

arguments (Input)
    spsd struct
    sensitivity_curve {mustBeMatrix}
    fill_value {mustBeNumeric}
end

arguments (Output)
    spsd_cal struct
end

% check if 'freq' exists in spsd
if ~isfield(spsd, "freq")
    error('"freq" is missing in spsd!')
end

spsd_cal = spsd;

% interpolate calibration
calibration = interp1(sensitivity_curve(:,1), ...
    sensitivity_curve(:,2), spsd.freq, "linear");
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
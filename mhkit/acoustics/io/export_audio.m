function export_audio(filename, pressure_data, peak_voltage, gain, resample_multiplier)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Export human-scaled audio data from hydrophone recording into WAV file
%
% Parameters
% ------------
%   filename: string
%       Output filename for the WAV file (without extension)
%   pressure_data: struct
%       Sound pressure or voltage data structure containing:
%       pressure_data.data : Pressure/Voltage data array
%       pressure_data.sensitivity : Sensitivity of the hydrophone [dB re 1 V/uPa]
%       pressure_data.fs : Sampling frequency [Hz]
%       pressure_data.units : Data units ('Pa' or 'V')
%       pressure_data.peak_voltage : (Optional) Peak voltage
%       pressure_data.valid_max : (Optional) Valid maximum voltage limit
%   peak_voltage: float (optional)
%       Peak voltage of the analog-to-digital converter. Default is empty.
%   gain: float (optional)
%       Gain to multiply the original time series by. Default is 1.
%   resample_multiplier: int (optional)
%       Multiplier for resampling the pressure to speed up the recording,
%       which is useful for listening to low frequency sound. Default is 1.
%
% Returns
% ---------
%   None (writes WAV file to disk)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    filename {mustBeText}
    pressure_data struct
    peak_voltage = []
    gain {mustBeNumeric} = 1
    resample_multiplier {mustBeInteger, mustBeGreaterThanOrEqual(resample_multiplier, 1)} = 1
end

% Validate structure fields
if ~isfield(pressure_data, 'data')
    error('MHKiT:acoustics:export_audio:MissingField', ...
        'pressure_data structure must contain data field');
end

if ~isfield(pressure_data, 'fs')
    error('MHKiT:acoustics:export_audio:MissingField', ...
        'pressure_data structure must contain fs field');
end

if ~isfield(pressure_data, 'sensitivity')
    error('MHKiT:acoustics:export_audio:MissingField', ...
        'pressure_data structure must contain sensitivity field');
end

fs = double(pressure_data.fs);

% Handle peak_voltage resolution
if isempty(peak_voltage)
    if isfield(pressure_data, 'peak_voltage')
        peak_voltage = double(pressure_data.peak_voltage);
    elseif isfield(pressure_data, 'valid_max')
        peak_voltage = double(pressure_data.valid_max);
    else
        error('MHKiT:acoustics:export_audio:MissingPeakVoltage', ...
            'pressure_data must contain a peak_voltage/valid_max field or peak_voltage must be supplied.');
    end
end

% Prepare signal data
v = double(pressure_data.data);

% Resample/Speed up recording if resample_multiplier > 1
if resample_multiplier > 1
    % Keep fs unchanged, but speed up by skipping or resampling data points.
    % To match Python's pandas.interp behavior:
    num_samples = length(v);
    new_num_samples = floor(num_samples / resample_multiplier);
    
    if new_num_samples < 2
        error('MHKiT:acoustics:export_audio:ResampleError', ...
            'Signal is too short for the specified resample_multiplier.');
    end
    
    % Interpolate the signal to a shorter coordinate series
    old_indices = 1:num_samples;
    new_indices = linspace(1, num_samples, new_num_samples);
    v = interp1(old_indices, v, new_indices, 'linear')';
end

% Convert from Pascals to Voltage if units are "Pa"
has_units = isfield(pressure_data, 'units');
if ~has_units || strcmpi(pressure_data.units, 'Pa')
    upa = v * 1e6;
    v = upa * 10^(double(pressure_data.sensitivity) / 20);
end

% Scale and Normalize
max_val = max(abs(v));
denom = min(max_val * double(gain), double(peak_voltage));
if denom > 0
    v = v / denom;
else
    warning('MHKiT:acoustics:export_audio:ZeroSignal', ...
        'Signal has zero amplitude or invalid normalization divisor.');
end

% Ensure output file has .wav extension
[filepath, name, ext] = fileparts(filename);
if isempty(ext)
    output_filename = fullfile(filepath, strcat(name, '.wav'));
else
    output_filename = filename;
end

% Write as 16-bit WAV file
try
    audiowrite(output_filename, v, fs, 'BitsPerSample', 16);
    fprintf('Audio exported successfully to: %s\n', output_filename);
catch ME
    error('MHKiT:acoustics:export_audio:WriteError', ...
        'Failed to write audio file: %s', ME.message);
end

end

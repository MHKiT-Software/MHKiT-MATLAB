function export_audio(filename, pressure_data, gain)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Export human-scaled audio data from hydrophone recording into WAV file
%
% Parameters
% ------------
%   filename: string
%       Output filename for the WAV file (without extension)
%   pressure_data: struct
%       Sound pressure data structure containing:
%       pressure_data.data : Pressure data array [Pa]
%       pressure_data.sensitivity : Sensitivity of the hydrophone [dB re 1 V/uPa]
%       pressure_data.fs : Sampling frequency [Hz]
%   gain: float
%       Gain to multiply the original time series by. Default: 1
%
% Returns
% ---------
%   None (writes WAV file to disk)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    filename {mustBeText}
    pressure_data struct
    gain {mustBeNumeric} = 1
end

% Validate pressure_data structure has required fields
if ~isfield(pressure_data, 'data')
    error('MHKiT:acoustics:export_audio:MissingField', ...
        'pressure_data structure must contain data field');
end

if ~isfield(pressure_data, 'sensitivity')
    error('MHKiT:acoustics:export_audio:MissingField', ...
        'pressure_data structure must contain sensitivity field');
end

if ~isfield(pressure_data, 'fs')
    error('MHKiT:acoustics:export_audio:MissingField', ...
        'pressure_data structure must contain fs field');
end

% Validate field types
if ~isnumeric(pressure_data.data)
    error('MHKiT:acoustics:export_audio:InvalidInput', ...
        'pressure_data.data must be numeric');
end

if ~isnumeric(pressure_data.sensitivity) || ~isscalar(pressure_data.sensitivity)
    error('MHKiT:acoustics:export_audio:InvalidInput', ...
        'pressure_data.sensitivity must be a numeric scalar');
end

if ~isnumeric(pressure_data.fs) || ~isscalar(pressure_data.fs) || pressure_data.fs <= 0
    error('MHKiT:acoustics:export_audio:InvalidInput', ...
        'pressure_data.fs must be a positive numeric scalar');
end

if ~isnumeric(gain) || ~isscalar(gain)
    error('MHKiT:acoustics:export_audio:InvalidInput', ...
        'gain must be a numeric scalar');
end

% Process audio data following Python implementation
% Convert from Pascals to microPascals
upa = pressure_data.data * 1e6;

% Convert to voltage waveform using sensitivity
v = upa * 10^(pressure_data.sensitivity / 20);  % in V

% Normalize to maximum absolute value and apply gain
max_abs_v = max(abs(v));
if max_abs_v > 0
    v = v / max_abs_v * gain;
else
    warning('MHKiT:acoustics:export_audio:ZeroSignal', ...
        'Signal has zero amplitude - no normalization applied');
end

% Ensure filename has .wav extension
[filepath, name, ext] = fileparts(filename);
if isempty(ext)
    output_filename = fullfile(filepath, strcat(name, '.wav'));
else
    output_filename = filename;
end

% Write as 16-bit WAV file
try
    audiowrite(output_filename, v, pressure_data.fs, 'BitsPerSample', 16);
    fprintf('Audio exported successfully to: %s\n', output_filename);
catch ME
    error('MHKiT:acoustics:export_audio:WriteError', ...
        'Failed to write audio file: %s', ME.message);
end

end

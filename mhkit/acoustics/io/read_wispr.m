function out = read_wispr(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read WISPR .dat file and return voltage time series structure.
%
% Parameters
% ------------
%   filename: string
%       Path to WISPR .dat file.
%
% Returns
% ---------
%   out: struct
%       Contains Sound pressure [Pa] or Voltage [V] along with metadata
%       out.data : Time series data [V]
%       out.time : Time vector
%       out.units : Data units ('V')
%       out.resolution : Voltage resolution [V]
%       out.valid_min : Minimum valid voltage [V]
%       out.valid_max : Maximum valid voltage [V]
%       out.fs : Sampling frequency [Hz]
%       out.filename : Original filename stem
%       out.gain : Gain setting in dB (6 dB intervals)
%       out.peak_voltage : Peak voltage [V]
%       out.file_length_sec : Total file length in seconds
%       out.instrument_id : Instrument ID of the WISPR sensor
%       out.sfw_version : WISPR software version
%       out.location_id : Location ID of deployment
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    filename {mustBeText}
end

arguments (Output)
    out struct
end

% Read the metadata
metadata = read_wispr_metadata(filename);

% Clean up metadata
start_time = datetime(metadata.time, 'InputFormat', 'MM:dd:yy:HH:mm:ss');
fs = double(metadata.sampling_rate);
peak_voltage = double(metadata.adc_vref);
bits_per_sample = double(metadata.sample_size * 8);

% Read binary data from wispr file
% Data is recorded in 16 or 24-bit by the ADC, saved in 32-bit format by the microcontroller,
% and finally converted to 16-bit within the WISPR code before being written to file.
if bits_per_sample == 24
    % 24-bit data is stored as 24-bit signed integers in little-endian format
    data = read_24bit_data(filename, true, '<');
else
    fid = fopen(filename, 'rb', 'ieee-le');
    if fid == -1
        error('MHKiT:acoustics:read_wispr:FileNotFound', 'Could not open file: %s', filename);
    end
    % skip header lines
    fseek(fid, 512, 'bof');
    data = fread(fid, Inf, 'int16=>double');
    fclose(fid);
end

% Normalize and then scale to peak voltage
max_count = 2^(bits_per_sample - 1);
% Use 64 bit float for decimal accuracy
raw_voltage = double(data) / max_count * peak_voltage;

% Establish time vector
dt = seconds(metadata.file_length_sec) / length(data);
time_vector = start_time + (0:double(length(data))-1)' * dt;

% Construct returned struct matching standard MHKiT acoustics structure
out = struct();
out.data = raw_voltage;
out.time = time_vector;
out.units = 'V';
out.resolution = round(peak_voltage / max_count, 9);
out.valid_min = -peak_voltage;
out.valid_max = peak_voltage;
out.fs = fs;

[~, name_stem, ~] = fileparts(filename);
out.filename = name_stem;

if isfield(metadata, 'gain')
    out.gain = metadata.gain * 6;
else
    out.gain = [];
end

out.peak_voltage = peak_voltage;

if isfield(metadata, 'file_length_sec')
    out.file_length_sec = metadata.file_length_sec;
else
    out.file_length_sec = [];
end

if isfield(metadata, 'instrument_id')
    out.instrument_id = metadata.instrument_id;
else
    out.instrument_id = [];
end

if isfield(metadata, 'version')
    out.sfw_version = metadata.version;
else
    out.sfw_version = [];
end

if isfield(metadata, 'location_id')
    out.location_id = metadata.location_id;
else
    out.location_id = [];
end

end

function data_int32 = read_24bit_data(filename, is_signed, endian)
% READ_24BIT_DATA Reads 24-bit data from a binary file into a double column.
%
%   filename : Path of the file
%   is_signed : logical flag indicating if signed conversion is required
%   endian : Character string ('<' for little-endian, '>' for big-endian)

if nargin < 2
    is_signed = true;
end
if nargin < 3
    endian = '<';
end

fid = fopen(filename, 'rb');
if fid == -1
    error('MHKiT:acoustics:read_wispr:FileNotFound', 'Could not open file: %s', filename);
end
raw_bytes = fread(fid, Inf, 'uint8=>uint8');
fclose(fid);

% Ensure the file size is a multiple of 3 bytes
if mod(length(raw_bytes), 3) ~= 0
    error('MHKiT:acoustics:read_wispr:InvalidSize', 'File size is not a multiple of 3 bytes (24 bits)');
end

% Reshape into 3 rows of bytes matching the 3 bytes per sample
data_3bytes = reshape(raw_bytes, 3, []);

if endian == "<"
    % Little-endian byte combination
    val = double(data_3bytes(1, :)) + ...
          double(data_3bytes(2, :)) * 256 + ...
          double(data_3bytes(3, :)) * 65536;
else
    % Big-endian byte combination
    val = double(data_3bytes(3, :)) + ...
          double(data_3bytes(2, :)) * 256 + ...
          double(data_3bytes(1, :)) * 65536;
end

if is_signed
    % Perform 24-bit to 32-bit sign extension mathematically
    neg_idx = val >= 8388608; % 2^23
    val(neg_idx) = val(neg_idx) - 16777216; % 2^24
end

data_int32 = val';

end

function out = read_iclisten(filename, sensitivity, use_metadata)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read .wav file from an icListen hydrophone
%
% Returns voltage timeseries if sensitivity not provided, returns pressure timeseries if it is provided
%
% Parameters
% ------------
%   filename: string
%       Input filename
%   sensitivity: float
%       Hydrophone calibration sensitivity in dB re 1 V/uPa
%       Should be negative. Default: None
%   use_metadata: logical
%       If true and 'sensitivity' is empty, applies sensitivity value
%       stored in the .wav file. If false and 'sensitivity' is empty, a
%       sensitivity value is not applied
%
% Returns
% ---------
%   out: struct
%       Contains Sound pressure [Pa] or Voltage [V] along with metadata
%       out.data : Time series data [Pa] or [V]
%       out.time : Time vector
%       out.units : Data units ('Pa' or 'V')
%       out.fs : Sampling frequency [Hz]
%       out.filename : Original filename
%       out.serial_number : Device serial number
%       out.bits_per_sample : Bits per sample
%       out.peak_voltage : Peak voltage [V]
%       out.humidity : Humidity reading
%       out.temperature : Temperature reading
%       out.accelerometer : Accelerometer data
%       out.magnetometer : Magnetometer data
%       out.count_at_peak_voltage : Count at peak voltage
%       out.sequence_num : Sequence number
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    filename {mustBeText}
    sensitivity = []
    use_metadata logical = true
end

arguments (Output)
    out struct
end

% decode metadata
metaraw = audioinfo(filename);
% get start_time
slen = strlength(metaraw.Title);
stime = datetime(metaraw.Title(slen-14:end),"InputFormat",'yyyyMMdd_HHmmss');
start_time = string(stime,'yyyy-MM-dd HH:mm:ss');

meta = split(metaraw.Comment, ', ');
% get peak voltage
temp = split(meta{1},' ');
peak_voltage = str2double(temp{1});
if isempty(sensitivity) & use_metadata
    % get sensitivity
    temp = split(meta{2},' ');
    sensitivity = str2double(temp{1});
end

out = read_hydrophone(filename,peak_voltage,sensitivity, 0, start_time);

% add extra metadata
out.serial_number = metaraw.Artist;
out.bits_per_sample = metaraw.BitsPerSample;
out.peak_voltage = peak_voltage;
out.humidity = meta{3};
out.temperature = meta{4};
out.accelerometer = meta{5};
out.magnetometer = meta{6};
out.count_at_peak_voltage = meta{7};
out.sequence_num = meta{8};

end

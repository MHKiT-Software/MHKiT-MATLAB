function out = read_soundtrap(filename, sensitivity, gain)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read .wav file from an Ocean Instruments SoundTrap hydrophone
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
%   gain: float
%       Amplifier gain in dB re 1 V/μPa
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
%       out.valid_min : Minimum valid data value
%       out.valid_max : Maximum valid data value
%       out.sensitivity : Sensitivity value (if provided)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    filename {mustBeText}
    sensitivity = []
    gain {mustBeNumeric} = 0
end

arguments (Output)
    out struct
end

% parse start time from filename
s1 = strsplit(filename, '.');
s1 = s1(end-1);
stime = datetime(s1, InputFormat="yyMMddHHmmss");
start_time = string(stime, 'yyyy-MM-dd HH:mm:ss');

peak_voltage = 1; % soundtrap uses peak voltage of 1 V

out = read_hydrophone(filename, peak_voltage, sensitivity, gain, start_time);

end

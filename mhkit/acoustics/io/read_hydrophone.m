function out = read_hydrophone(filename, peak_voltage, sensitivity, gain, start_time)
%{
    Read .wav file from a hydrophone. Returns voltage timeseries if sensitivity not
    provided, returns pressure timeseries if it is provided.

    Parameters
    ----------
    filename: string
        Input filename
    peak_voltage: double
        Peak voltage supplied to the analog to digital converter (ADC) in V.
        (Or 1/2 of the peak to peak voltage).
    sensitivity: double
        Hydrophone calibration sensitivity in dB re 1 V/uPa.
        Should be negative. Default: None.
    gain: double
        Amplifier gain in dB re 1 V/uPa. Default 0.
    start_time: string
        Start time in the format yyyy-MM-ddTHH:mm:ss

    Returns
    -------
    out: struct
        Contains Sound pressure [Pa] or Voltage [V] along with metadata
%}

arguments (Input)
    filename {mustBeText}
    peak_voltage {mustBeNumeric}
    sensitivity = []
    gain {mustBeNumeric} = 0
    start_time {mustBeText} = "2024-01-01T00:00:00"
end

arguments (Output)
    out struct
end

[y, fs] = audioread(filename); % read in file
v = y * peak_voltage; % scale based on peak voltage

out = struct();
len = length(v) / fs;
start = datetime(start_time);
ending = start + seconds(len);
out.time = linspace(start, ending, length(v));

if ~isempty(sensitivity)
    if sensitivity > 0
        error('Sensitivity must be negative numeric input!')
    end
    sense = sensitivity - gain;
    Sf = 10^(sense/20);
    Pu = v / Sf;  % sound pressure in uPa
    P = Pu/10^6;  % sound pressure in Pa
    out.data = P;
    out.units = 'Pa';
    out.sensitivity = Sf;
    out.valid_min = -peak_voltage / Sf / 1e6;
    out.valid_max = peak_voltage / Sf / 1e6;
else
    out.data = v; % scaled voltage output
    out.units = 'V';
    out.valid_min = -peak_voltage;
    out.valid_max = peak_voltage;
end

% add metadata
out.fs = fs;
[~,ff,ext] = fileparts(filename);
out.filename = strcat(ff,ext);

end
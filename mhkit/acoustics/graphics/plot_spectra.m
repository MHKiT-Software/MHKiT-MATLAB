function plot_spectra(spsdl, fmin, fmax)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot spectra of the sound pressure spectral density level
%
% Parameters
% ------------
%   spsdl: struct
%       Mean square sound pressure spectral density level in dB rel 1 uPa^2/Hz
%       spsdl.data : Spectral density level data [dB rel 1 uPa^2/Hz]
%       spsdl.freq : Frequency vector [Hz]
%       spsdl.time : Time vector
%       spsdl.name : Data name
%       spsdl.units : Data units
%   fmin: float
%       Lower frequency band limit (lower limit of the hydrophone) [Hz]
%   fmax: float
%       Upper frequency band limit (Nyquist frequency) [Hz]
%
% Returns
% ---------
%   This function creates a figure but does not return a value
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    spsdl struct
    fmin {mustBeNumeric} = 10
    fmax {mustBeNumeric} = 100000
end

% Validate spsdl structure
validate_spsdl_struct(spsdl, 'plot_spectra', ...
    'required_fields', {{'data', 'freq', 'name', 'units'}}, ...
    'check_dimensions', false);

% check fmax
fn = max(spsdl.freq);
fmax = fmax_warning(fn, fmax);

% select freq range
idxf = find(spsdl.freq>=fmin & spsdl.freq<=fmax);

figure;
plot(spsdl.freq(idxf), spsdl.data(idxf,:))
xscale log;
xlim([fmin fmax])
xlabel('Frequency [Hz]'); 
ylabel({spsdl.name;strcat('[',spsdl.units,']')});
grid on;

end

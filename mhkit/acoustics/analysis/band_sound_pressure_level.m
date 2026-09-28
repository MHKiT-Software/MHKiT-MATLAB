function mspl = band_sound_pressure_level(spsd, octave, base, fmin, fmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate band-averaged sound pressure levels from the mean square sound pressure spectral density (SPSD)
%
% Parameters
% ------------
%   spsd: struct
%       Mean square sound pressure spectral density in [Pa^2/Hz]
%       spsd.data : Spectral density data [Pa^2/Hz]
%       spsd.freq : Frequency vector [Hz]
%       spsd.time : Time vector
%       spsd.fs : Sampling frequency [Hz]
%   octave: float
%       Octave subdivision (1 = full octave, 3 = third-octave, etc.)
%   base: float
%       Octave base subdivision (2 = true octave, 10 = decade octave, etc.)
%   fmin: float
%       Lower frequency band limit (lower limit of the hydrophone) [Hz]
%   fmax: float
%       Upper frequency band limit (Nyquist frequency) [Hz]
%
% Returns
% ---------
%   mspl: struct
%       Sound pressure level [dB re 1 uPa] indexed by time and frequency of specified bandwidth
%       mspl.data : Sound pressure level data [dB re 1 uPa]
%       mspl.freq : Center frequencies of bands [Hz]
%       mspl.time : Time vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    spsd struct
    octave {mustBeInteger}
    base {mustBeInteger} = 2
    fmin {mustBeNumeric} = 10
    fmax {mustBeNumeric} = 100000
end

arguments (Output)
    mspl struct
end

% Validate spsd structure
validate_spsd_struct(spsd, 'band_sound_pressure_level', ...
    'required_fields', {{'freq', 'time', 'fs'}});

reference = 1e-12; % Pa^2, = 1 uPa^2

fmax = fmax_warning(spsd.fs/2, fmax);

% Create frequency bands
[obins, band] = create_frequency_bands(octave, base, fmin, fmax);

nBands = numel(band.center_freq);
nTime = numel(spsd.time);
mspl.data = zeros(nBands, nTime);

for i = 1:nBands
    band_range = [band.lower_limit(i), band.upper_limit(i)];
    idx = find(spsd.freq >= band_range(1) & spsd.freq <= band_range(2));
    if numel(idx) < 2
        % Interpolate if only one frequency in band
        spsd_slc = interp1(spsd.freq, spsd.data, band_range, 'linear', 'extrap');
        x = band_range;
    else
        spsd_slc = spsd.data(idx, :);
        x = spsd.freq(idx);
    end
    % Integrate spectral density by frequency for each time
    for t = 1:nTime
        mspl.data(i, t) = 10 * log10(trapz(x, spsd_slc(:, t)) / reference);
    end
end

mspl.freq = band.center_freq;
mspl.time = spsd.time;

end

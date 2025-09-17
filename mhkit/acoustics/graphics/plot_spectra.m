function plot_spectra(spsdl, fmin, fmax)
%{
    Plots the spectrogram of the sound pressure spectral density level.

    Parameters
    ----------
    spsdl: struct
        Mean square sound pressure spectral density level in dB rel 1 uPa^2/Hz
    fmin: integer
        Lower frequency band limit (lower limit of the hydrophone). Default: 10 Hz
    fmax: integer
        Upper frequency band limit (Nyquist frequency). Default: 100000 Hz

    Returns
    -------
    sp: figure
%}

arguments (Input)
    spsdl struct
    fmin {mustBeNumeric} = 10
    fmax {mustBeNumeric} = 100000
end

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
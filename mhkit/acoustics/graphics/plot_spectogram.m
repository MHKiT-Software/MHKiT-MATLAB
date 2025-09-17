function plot_spectogram(spsdl, fmin, fmax, cm, vmin, vmax)
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
    cm {mustBeText} = 'hot'
    vmin {mustBeNumeric} = 0
    vmax {mustBeNumeric} = 0
end

% check fmax
fn = max(spsdl.freq);
fmax = fmax_warning(fn, fmax);

% select freq range
idxf = find(spsdl.freq>=fmin & spsdl.freq<=fmax);

figure;
pcolor(spsdl.time, spsdl.freq(idxf), spsdl.data(idxf,:))
shading interp
yscale log;
colormap(cm);
cb = colorbar();
if ~isequal(vmax, 0)
    clim([vmin vmax]);
end
xlabel('Time'); ylabel('Frequency [Hz]')
cb.Label.String = spsdl.units;

end
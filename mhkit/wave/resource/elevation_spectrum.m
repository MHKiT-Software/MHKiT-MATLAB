function S=elevation_spectrum(ts,sample_rate,nnft,time,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates wave spectra from wave probe timeseries
%
% Parameters
% ------------
%     ts: matrix or table
%           Wave probe time-series data, with each column a different time
%           series
%
%     sample_rate: float
%           Data frequency (Hz)
%
%     nnft: integer
%           Number of bins in the Fast Fourier Transform
%
%     time: vector or table
%           time (s)
%
%     window: string (Optional)
%        Signal window type. "hamming" is used by default given the broadband
%        nature of waves. See scipy.signal.get_window for more options.
%        to call: elevation_spectrum(ts,sample_rate,nnft,time,"window",window)
%
%     detrend: logical (Optional)
%        Specifies if a linear trend is removed from the data before calculating
%        the wave energy spectrum.  Data is detrended by default.
%        to call: elevation_spectrum(ts,sample_rate,nnft,time,"detrend",detrend)
%
%     noverlap: integer (Optional)
%       Number of points to overlap between segments. If None,
%       noverlap = nperseg / 2.  Defaults to None.
%       to call: elevation_spectrum(ts,sample_rate,nnft,time,"noverlap",noverlap)
%
% Returns
% ---------
%     S: structure
%
%
%         S.spectrum: vector or matrix Spectral Density (m^2/Hz) per probe
%
%         S.type: 'Spectra from Time Series'
%
%         S.frequency: frequency [Hz]
%
%         S.sample_rate: sample_rate
%
%         S.nnft: nnft
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    ts
    sample_rate
    nnft
    time
    options.window = "hamming";
    options.detrend = true;
    options.noverlap = py.None;
end

ts = py.pandas.DataFrame(py.numpy.array(ts));

nnft = int32(nnft);

spectra_py = py.mhkit.wave.resource.elevation_spectrum( ...
    ts, sample_rate, nnft, ...
    pyargs( ...
        'window', options.window, ...
        'detrend', options.detrend, ...
        'noverlap', options.noverlap ...
    ) ...
);

spectra_py = typecast_from_mhkit_python(spectra_py);

S = struct();

S.type = 'Spectra from Timeseries';
S.spectrum = spectra_py.data;
S.frequency = spectra_py.index.data;
S.nnft = nnft;
S.sample_rate = sample_rate;
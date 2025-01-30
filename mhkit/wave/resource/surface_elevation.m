function wave_elevation=surface_elevation(S,time_index,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates wave elevation time series from spectrum using a random phase
%
% Parameters
% ------------
%    S: Spectral Density (m^2/Hz)
%       Pandas data frame
%           To make a pandas data frame from user supplied frequency and spectra
%           use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(frequency,spectra)
%
%       OR
%
%       structure of form:
%           S.spectrum: Spectral Density (m^2/Hz)
%
%           S.type: String of the spectra type, i.e. Bretschneider,
%           time series, date stamp etc.
%
%           S.frequency: frequency (Hz)
%
%    time_index: array
%        Time used to create the wave elevation time series [s]
%
%    seed: Int (optional)
%        random seed
%        to call: wave_elevation(S,time_index,"seed",seed)
%
%    frequency_bins: vector (optional)
%       Bin widths for frequency of S. Required for unevenly sized bins
%       to call: wave_elevation(S,time_index,"frequency_bins",frequency_bins)
%
%    phases: vector or matrix (optional)
%       Explicit phases for frequency components (overrides seed)
%       to call: wave_elevation(S,time_index,"phases",phases)
%
%    method: str (optional)
%       Method used to calculate the surface elevation. 'ifft' (Inverse Fast Fourier
%       Transform) used by default if the given frequency_bins==None. 'sum_of_sines'
%       explicitly sums each frequency component and used by default if
%       frequency_bins are provided. The 'ifft' method is significantly faster.
%
% Returns
% ---------
%    wave_elevation: structure
%
%         wave_elevation.elevation: Wave surface elevation (m)
%
%         wave_elevation.type: 'Time Series from Spectra'
%
%         wave_elevation.time
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    S
    time_index
    options.seed {mustBeNumeric} = 123;
    options.frequency_bins = py.None;
    options.phases = py.None;
    options.method = "ifft";
end

S_py = typecast_spectra_to_mhkit_python(S);

frequency = S.frequency ;

if (isa(time_index,'py.numpy.ndarray') ~= 1)

    time_index = py.numpy.array(time_index);

end

if (isa(options.frequency_bins,'py.NoneType')~=1)
    if isnumeric(options.frequency_bins)

        options.frequency_bins = py.numpy.array(options.frequency_bins);
    else
        ME = MException('MATLAB:significant_wave_height','frequency_bins need to be of numeric type');
        throw(ME);
    end
end

if (isa(options.phases,'py.NoneType')~=1)
    if isnumeric(options.phases)
        options.phases = py.numpy.array(options.phases);
    else
        ME = MException('MATLAB:significant_wave_height','phases need to be of numeric type');
        throw(ME);
    end

end

if ~(strcmp(options.method, "ifft") || strcmp(options.method, "sum_of_sines"))
    ME = MException('MATLAB:significant_wave_height','Invalid method. method should be either "ifft" or "sum_of_sines"');
    throw(ME);
end

eta_py=py.mhkit.wave.resource.surface_elevation(S_py, time_index);

eta = typecast_from_mhkit_python(eta_py);

wave_elevation.type='Time Series from Spectra';
wave_elevation.time=eta.index.data;
wave_elevation.elevation = eta.data;

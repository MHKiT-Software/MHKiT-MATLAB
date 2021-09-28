function l=wave_length(k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates wave length from wave number
% To compute: 2*pi/wavenumber
%    
% Parameters
% ------------
%    k: wave number (1/m)
%        structure of form:
%           k.values= wave number
%
%           k.frequency= frequency (Hz)
%
% Returns
% -------
%    l: double or array
%         Wave length [m] indexed by frequency
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

py.importlib.import_module('mhkit');
py.importlib.import_module('mhkit_python_utils');

if (isa(k,'py.pandas.core.frame.DataFrame')~=1)
    if (isstruct(k)==1)
        k=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(k.frequency,py.numpy.array(k.values),1);

    else
        ME = MException('MATLAB:wave_length','k needs to be a structure or Pandas dataframe, use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas to create one');
        throw(ME);
    end
end

l=double(py.array.array('d',py.numpy.nditer(py.mhkit.wave.resource.wave_length(k))));
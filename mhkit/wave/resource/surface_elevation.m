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
% Returns
% ---------
%    wave_elevation: structure
%
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
end

py.importlib.import_module('mhkit');
% py.importlib.import_module('numpy');
py.importlib.import_module('mhkit_python_utils');

frequency= S.frequency ;

if (isa(time_index,'py.numpy.ndarray') ~= 1)

    time_index = py.numpy.array(time_index);

end

if (isa(S,'py.pandas.core.frame.DataFrame')~=1)
    if (isstruct(S)==1)
        x=size(S.spectrum);
        li=py.list();
        if x(2)>1
            for i = 1:x(2)
                app=py.list(S.spectrum(:,i));
                li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);

            end
            S=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(frequency(:,1),li,int32(x(2)));
        elseif x(2)==1
            S=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(frequency,py.numpy.array(S.spectrum),int32(x(2)));
        end
    else
        ME = MException('MATLAB:significant_wave_height','S needs to be a structure or Pandas dataframe, use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas to create one');
        throw(ME);
    end
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
     x=size(options.phases);
     li=py.list();
     if x(2)>1
         for i = 1:x(2)
             app=py.list(options.phases(:,i));
             li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);

         end
         options.phases=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(frequency(:,1),li,int32(x(2)));
      elseif x(2)==1
         options.phases=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(frequency,py.numpy.array(options.phases),int32(x(2)));
     end

    else
        ME = MException('MATLAB:significant_wave_height','phases need to be of numeric type');
        throw(ME);
    end

end


eta=py.mhkit.wave.resource.surface_elevation(S,time_index,pyargs('seed',...
    py.int(options.seed),'frequency_bins',options.frequency_bins,'phases',options.phases));


vals=double(py.array.array('d',py.numpy.nditer(eta.values)));
 sha=cell(eta.values.shape);
 x=int64(sha{1,1});
 y=int64(sha{1,2});
 vals=reshape(vals,[x,y]);

si=size(vals);

wave_elevation.elevation=vals;
% for i=1:si(2)
%    wave_elevation.spectrum{i}=vals(:,i);
% end

wave_elevation.type='Time Series from Spectra';

wave_elevation.time=double(py.array.array('d',py.numpy.nditer(eta.index)));


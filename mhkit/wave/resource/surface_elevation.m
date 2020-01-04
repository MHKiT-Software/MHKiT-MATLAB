function wave_elevation=surface_elevation(S,time_index,varargin)

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
%        Time used to create the wave elevation time series [s],
%        
%    seed: Int (optional)
%        random seed
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


py.importlib.import_module('mhkit');
py.importlib.import_module('numpy');
py.importlib.import_module('mhkit_python_utils');

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
            S=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(S.frequency(:,1),li,int32(x(2)));
        elseif x(2)==1
            S=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(S.frequency,py.numpy.array(S.spectrum),int32(x(2)));
        end
    else
        ME = MException('MATLAB:significant_wave_height','S needs to be a structure or Pandas dataframe, use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas to create one');
        throw(ME);
    end
end

 if nargin == 4 
     seed=varagin{1};
     eta=py.mhkit_wave_resource.surface_elevation(S,time_index,pyargs('seed',seed));
 else
     eta=py.mhkit.wave.resource.surface_elevation(S,time_index);
 end

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

    
    

function Cg=wave_celerity(k,h,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates wave celerity (group velocity)
%    
% Parameters
% ------------
%    k: wave number (1/m)
%       Pandas data frame
%           To make a pandas data frame from user supplied frequency and spectra
%           use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(frequency,spectra)
%
%        OR
%
%        structure of form:
%           k.values= wave number
%
%           k.frequency= frequency (Hz)
%
%    h: float
%         Water depth (m)
%
%    g: float (optional)
%         gravitational acceleration (m/s^2)
%         
%
% Returns
% -------
%     Cg: structure
% 
%
%       Cg.values: water celerity
%
%       Cg.frequency [Hz]
%
%       Cg.h: height [m]
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


py.importlib.import_module('mhkit');
py.importlib.import_module('mhkit_python_utils');

if (isa(k,'py.pandas.core.frame.DataFrame')~=1)
    if (isstruct(k)==1)
%         x=size(k.values);
%         li=py.list();
%         if x(2)>1 
%             for i = 1:x(2)
%                 app=py.list(k.values(:,i));
%                 li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);
%             
%             end
%             k=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(k.frequency,li,int32(x(2)));
%         elseif x(2)==1
            k=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(k.frequency,py.numpy.array(k.values),1);
        %end
    else
        ME = MException('MATLAB:wave_celerity','k needs to be a structure or Pandas dataframe, use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas to create one');
        throw(ME);
    end
end

if nargin == 3 
    Cgdf=py.mhkit.wave.resource.wave_celerity(k,h,pyargs('g',varargin{1}));
elseif nargin == 2
    Cgdf=py.mhkit.wave.resource.wave_celerity(k,h);
else
    ME = MException('MATLAB:wave_celerity','incorrect number of arguments');
        throw(ME);
end


Cg.values=double(py.array.array('d',py.numpy.nditer(Cgdf.values)));
Cg.frequency=double(py.array.array('d',py.numpy.nditer(Cgdf.index)));
Cg.h=h;

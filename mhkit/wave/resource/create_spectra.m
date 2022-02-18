function S1=create_spectra(spectra_type,frequency,Tp,Hs,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates wave spectra of user specified type
%    
% Parameters
% ------------
%     spectraType: String of spectra type
%         Options are: 'pierson_moskowitz_spectrum',
%         'bretschneider_spectrum', or 'jonswap_spectrum'
%
%     Frequency: float
%         Wave frequency (Hz)
%
%     Tp: float
%         Peak Period (s)
%
%     Hs: float
%         Significant Wave Height (s)
%
%     gamma: float (optional)
%         only an optional parameter for 'jonswap_spectrum'
%     
% Returns
% ---------
%     S: structure
%
%
%         S.spectrum=Spectral Density (m^2/Hz)
%
%         S.type=String of the spectra type, i.e. Bretschneider, 
%         time series, date stamp etc. 
%
%         S.frequency= frequency (Hz)
%
%         S.Tp
%
%         S.Hs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


py.importlib.import_module('mhkit');
py.importlib.import_module('numpy');

if (isa(frequency,'py.numpy.ndarray') ~= 1)
    frequency = py.numpy.array(frequency);
end

if strcmp(spectra_type,'pierson_moskowitz_spectrum')
    if nargin ==4
        S=py.mhkit.wave.resource.pierson_moskowitz_spectrum(frequency,Tp,Hs);
    else
        ME = MException('MATLAB:create_spectra','incorrect number of arguments');
         throw(ME);
    end
    
elseif strcmp(spectra_type,'bretschneider_spectrum')
    if nargin ~= 4
         ME = MException('MATLAB:create_spectra','Hs is needed for bretschneider_spectrum');
         throw(ME);
    end
    S=py.mhkit.wave.resource.bretschneider_spectrum(frequency,Tp,Hs);

elseif strcmp(spectra_type,'jonswap_spectrum')
    if nargin < 4
         ME = MException('MATLAB:create_spectra','Hs is needed for jonswap_spectrum');
         throw(ME);
    elseif nargin == 4
        S=py.mhkit.wave.resource.jonswap_spectrum(frequency,Tp,Hs);
    elseif nargin == 5 
        S=py.mhkit.wave.resource.jonswap_spectrum(frequency,Tp,Hs,pyargs('gamma',varargin{1}));
    else
        ME = MException('MATLAB:create_spectra','to many input arguments');
         throw(ME);
    end

else 
    ME = MException('MATLAB:create_spectra','Not a Valid Spectrum Type');
    throw(ME);
end

S1.spectrum=double(py.array.array('d',py.numpy.nditer(S.values))).';
char_arr=char(S.index.values);
S1.type=spectra_type;
S1.frequency=double(py.array.array('d',py.numpy.nditer(S.index))).';
S1.Tp=Tp;
if nargin == 4 
    S1.Hs=Hs;
end




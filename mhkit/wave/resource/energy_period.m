function Te=energy_period(S,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Parameters
% ------------
%   S: Spectral Density (m^2/Hz)
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
%     frequency_bins: vector (optional)
%       Bin widths for frequency of S. Required for unevenly sized bins
%
%
% Returns
% ---------
%    Te: float
%        Wave energy Period (s)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% py.importlib.import_module('numpy');
py.importlib.import_module('mhkit');
py.importlib.import_module('mhkit_python_utils');

% assign feq_bin
if nargin == 2
    freq_bins = py.numpy.array(varargin{1});
elseif nargin == 1
    freq_bins = py.None;
else
    ME = MException('MATLAB:energy_period','Incorrect number of input arguments');
    throw(ME);
end

S = typecast_spectra_to_mhkit_python(S);

% % convert matlab structure to pandas.dataframe
% if (isa(S,'py.pandas.core.frame.DataFrame')~=1)
% 
%     if (isstruct(S)==1)
%             x=size(S.spectrum);
%             S_test = S;
%             li=py.list();
%             if x(2)>1
%                 for i = 1:x(2)
%                     app=py.list(double(S.spectrum(:,i)));
%                     li=py.mhkit_python_utils.pandas_dataframe.lis(li,app);
% 
%                 end
%                 % S=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas_v2(double(S.frequency(:,1)),S.spectrum(:,1),x(2));
%                 % S=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(double(S.frequency(:,1)),li,x(2));
%                 S = typecast_spectra_to_mhkit_python(double(S.frequency(:,1)));
% 
%                 %%
%             elseif x(2)==1
%                 % S=py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas(double(S.frequency),py.numpy.array(S.spectrum),x(2));
% 
%             end
% 
%     else
%         ME = MException('MATLAB:energy_period','S needs to be a structure or Pandas dataframe, use py.mhkit_python_utils.pandas_dataframe.spectra_to_pandas to create one');
%         throw(ME);
%     end
% end

% calculate enegery_period here
Te=py.mhkit.wave.resource.energy_period(S,pyargs('frequency_bins',freq_bins));
Te=double(Te.values);

end
function df = typecast_spectra_to_mhkit_python(S)
% Initialize empty python list for column names
column_names = py.list();

% Check and correct frequency orientation if needed
if size(S.frequency, 1) ~= 1
    warning('Input frequency is not row-wise (1×N). Current shape is (%d×%d). Converting to row vector.', ...
        size(S.frequency, 1), size(S.frequency, 2));
    S.frequency = S.frequency(:)';  % Convert frequency to row vector
end

% Now we know frequency is row-wise, so only reshape spectrum if needed
if size(S.spectrum, 1) > 1  % spectrum is column vector
    S.spectrum = S.spectrum(:)';  % Convert to row vector
end

% Get dimensions of spectrum matrix
[rows, cols] = size(S.spectrum);

% Create the frequency index
freq_index = py.numpy.array(S.frequency);

% Handle single or multiple spectra
if rows == 1
    % Single spectrum case
    spectrum_data = py.numpy.array(S.spectrum);
    % Use type as column name, default to 'spectrum' if not provided
    if isfield(S, 'type')
        column_names.append(py.str(S.type));
    else
        column_names.append(py.str('spectrum'));
    end
else
    % Multiple spectra case
    % Convert all columns to numpy arrays and store column names
    spectrum_data = py.numpy.zeros([rows, cols]);
    for i = 1:rows
        spectrum_data(:,i) = py.numpy.array(S.spectrum(:,i));
        % Generate column name for each spectrum
        if isfield(S, 'type') && length(S.type) >= i
            column_names.append(py.str(S.type{i}));
        else
            column_names.append(py.str(sprintf('spectrum_%d', i)));
        end
    end
end

% Create pandas DataFrame with frequency as index
df = py.pandas.DataFrame(data=spectrum_data, ...
    index=freq_index, ...
    columns=column_names);
end
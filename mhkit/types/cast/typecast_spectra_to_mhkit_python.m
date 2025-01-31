function df = typecast_spectra_to_mhkit_python(S)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Typecast input spectra struct into a MHKiT-Python compatible
% pandas DataFrame.
%
% Parameters
% ----------
% S : struct
%     A structure containing:
%     - frequency : (1D array) Frequency values (Hz)
%     - spectrum  : (1D or 2D matrix) Spectral density values (m^2/Hz)
%     - type (optional) : (cell array of strings) Names for each spectrum column,
%                        e.g., 'JONSWAP', 'Pierson Moskowitz', etc.
%
% Returns
% -------
% df : Python pandas.DataFrame
%     A Pandas DataFrame with frequencies (Hz) as the index and
%     spectral density values (m^2/Hz) as columns. Column names are either
%     provided by S.type or automatically generated as 'spectrum_N'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize empty Python list for column names
column_names = py.list();

% Validate that frequency is a 1D array and convert it to a row vector if needed
freq_dims = size(S.frequency);
if length(freq_dims) > 2 || (freq_dims(1) ~= 1 && freq_dims(2) ~= 1)
    error('Frequency must be a 1D array. Current dimensions are %s', mat2str(freq_dims));
end
freq_length = length(S.frequency);
if freq_dims(1) ~= 1
    S.frequency = S.frequency(:)';  % Ensure frequency is row-wise
end

% Get initial spectrum dimensions
[spec_rows, spec_cols] = size(S.spectrum);

% Handle different spectrum shapes and orientations
if spec_rows == spec_cols && spec_rows == freq_length
    % Case 1: Square matrix where rows = cols = frequency length
    warning('Spectrum is square (%d×%d) and matches frequency length. Assuming input orientation is correct.', ...
            spec_rows, spec_cols);
elseif spec_cols == freq_length
    % Case 2: Spectrum already correctly shaped (rows = M, cols = frequency length)
    % No transformation needed
elseif spec_rows == freq_length
    % Case 3: Spectrum has incorrect orientation (transpose needed)
    S.spectrum = S.spectrum';
    [spec_rows, spec_cols] = size(S.spectrum);
else
    % Case 4: Invalid dimensions that cannot be reconciled
    error('Spectrum dimensions (%d×%d) do not match frequency length (%d). Expected 1×%d or M×%d.', ...
          spec_rows, spec_cols, freq_length, freq_length, freq_length);
end

% Convert frequency array to Python NumPy array for Pandas index
freq_index = py.numpy.array(S.frequency);

% Convert spectrum data to NumPy array
if spec_rows == 1
    % Case: Single spectrum
    spectrum_data = py.numpy.array(S.spectrum);
    if isfield(S, 'type')
        column_names.append(py.str(S.type));  % Use provided type as column name
    else
        column_names.append(py.str('spectrum'));  % Default column name
    end
else
    % Case: Multiple spectra
    spectrum_data = py.numpy.zeros([spec_rows, spec_cols]);
    for i = 1:spec_rows
        spectrum_data(:, i) = py.numpy.array(S.spectrum(i, :));  % Convert each row
        if isfield(S, 'type') && length(S.type) >= i
            column_names.append(py.str(S.type{i}));  % Use provided type name
        else
            column_names.append(py.str(sprintf('spectrum_%d', i)));  % Generate name
        end
    end
end

% Create a Pandas DataFrame with frequency as the index
df = py.pandas.DataFrame(data=spectrum_data, ...
                         index=freq_index, ...
                         columns=column_names);

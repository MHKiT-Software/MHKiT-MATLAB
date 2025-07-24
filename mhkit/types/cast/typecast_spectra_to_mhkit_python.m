function df = typecast_spectra_to_mhkit_python(S)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Typecast input spectra struct into a MHKiT-Python compatible
% pandas DataFrame. This function handles conversion of MATLAB spectral data
% to a format compatible with MHKiT-Python's wave analysis functions.
%
% Parameters
% ----------
% S : struct
%     A structure containing:
%     - frequency : (1D array) Frequency values (Hz)
%     - spectrum  : (1D or 2D matrix) Spectral density values (m^2/Hz)
%                  For 2D matrices, each ROW corresponds to a different spectrum,
%                  and each COLUMN corresponds to a frequency point.
%
% Returns
% -------
% df : Python pandas.DataFrame
%     A Pandas DataFrame with:
%     - Index: frequency values (Hz)
%     - Data: spectral density values (m^2/Hz), where each row corresponds to
%            a frequency point and each column represents a different spectrum
%
% Example
% -------
% S.frequency = [0.01, 0.02, 0.03];  % Hz
% S.spectrum = [0.1, 0.2, 0.3;       % First spectrum
%               0.4, 0.5, 0.6];      % Second spectrum
% df = typecast_spectra_to_mhkit_python(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Validate and ensure frequency is a row vector for consistent handling
freq_dims = size(S.frequency);
if length(freq_dims) > 2 || (freq_dims(1) ~= 1 && freq_dims(2) ~= 1)
    error('Frequency must be a 1D array. Current dimensions are %s', mat2str(freq_dims));
end
freq_length = length(S.frequency);
if freq_dims(1) ~= 1
    S.frequency = S.frequency(:)';  % Convert to row vector if it's a column
end

% Get spectrum dimensions for orientation check
[spec_rows, spec_cols] = size(S.spectrum);

% Handle spectrum orientation
% The input spectrum should have rows matching the frequency length
% because each row represents a different spectrum
if spec_rows == freq_length
    % Spectrum is correctly oriented (frequencies in rows)
    spectrum_formatted = S.spectrum;
elseif spec_cols == freq_length
    % Transpose needed to get frequencies in rows
    spectrum_formatted = S.spectrum';
else
    error('Spectrum dimensions (%d×%d) do not match frequency length (%d).', ...
          spec_rows, spec_cols, freq_length);
end

% Convert MATLAB arrays to Python objects
freq_index = py.numpy.array(S.frequency);
spectrum_data = py.numpy.array(spectrum_formatted);

% Create DataFrame with frequency as index
% The resulting DataFrame will have frequencies as the index
% and each column will contain a different spectrum
df = py.pandas.DataFrame(data=spectrum_data, index=freq_index);
end

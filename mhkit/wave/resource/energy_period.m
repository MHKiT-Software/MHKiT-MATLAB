function Te = energy_period(S, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates wave energy period Te from wave spectra
%
% Parameters
% ------------
%   S: structure or numeric array
%       If structure:
%           S.spectrum   - Spectral Density (m^2/Hz)
%           S.frequency  - Frequency (Hz)
%       If numeric:
%           S is spectral density array (vector or matrix)
%           varargin{1} must contain frequency vector
%
%   frequency_bins: vector (optional)
%       Frequency bin widths [Hz]. Required for unevenly spaced bins.
%
% Returns
% ---------
%   Te: double
%       Wave energy period [s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Extract spectrum and frequency
    if isstruct(S)
        spectrum = S.spectrum;
        frequency = S.frequency;
    elseif isnumeric(S)
        if nargin < 2
            error('When S is numeric, frequency vector must be provided as second argument');
        end
        spectrum = S;
        frequency = varargin{1};
        varargin(1) = [];
    else
        error('Input S must be a struct or numeric array');
    end

    % Determine frequency bin widths
    if ~isempty(varargin)
        freq_bins = varargin{1};
    else
        df = diff(frequency);
        freq_bins = mean(df);
    end

    % Ensure proper shapes
    frequency = frequency(:);
    if isvector(spectrum)
        spectrum = spectrum(:);
    end

    % Check size consistency
    if length(frequency) ~= size(spectrum,1)
        error('Length of frequency vector must match number of rows in spectrum');
    end

    % Calculate moments: m0 and m-1
    if isscalar(freq_bins)
        m0 = sum(spectrum .* freq_bins, 1);
        m_neg1 = sum((spectrum ./ frequency) .* freq_bins, 1);
    else
        freq_bins = freq_bins(:);
        if length(freq_bins) ~= length(frequency)
            error('Length of freq_bins must match frequency vector');
        end
        m0 = sum(spectrum .* freq_bins, 1);
        m_neg1 = sum((spectrum ./ frequency) .* freq_bins, 1);
    end

    % Calculate energy period
    Te = m_neg1 ./ m0;

end

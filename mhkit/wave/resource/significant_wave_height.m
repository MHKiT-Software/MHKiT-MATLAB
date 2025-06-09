function Hm0 = significant_wave_height(S, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculates significant wave height Hm0 from spectra
%
% Parameters
% ------------
%     S: structure or numeric array
%         If structure:
%             S.spectrum: Spectral Density (m^2/Hz)
%             S.frequency: frequency (Hz)
%         If numeric:
%             S is assumed to be spectral density vector/matrix
%             varargin{1} must contain frequency vector
%
%     frequency_bins: vector (optional)
%         Bin widths for frequency of S. Required for unevenly sized bins
%
% Returns
% ---------
%     Hm0: double
%         Significant Wave Height (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        error('Input S must be either a struct with fields .spectrum and .frequency, or a numeric array');
    end

    % Extract frequency bin widths
    if ~isempty(varargin)
        freq_bins = varargin{1};
    else
        df = diff(frequency);
        freq_bins = mean(df);  % assume uniform bin width
    end

    % Ensure column vector format
    frequency = frequency(:);
    if isvector(spectrum)
        spectrum = spectrum(:);
    end

    % Check that frequency matches spectrum
    if length(frequency) ~= size(spectrum,1)
        error('Length of frequency vector must match number of rows in spectrum');
    end

    % Calculate zeroth moment m0
    if isscalar(freq_bins)
        m0 = sum(spectrum .* freq_bins, 1);
    else
        freq_bins = freq_bins(:);
        if length(freq_bins) ~= length(frequency)
            error('Length of freq_bins must match frequency vector');
        end
        m0 = sum(spectrum .* freq_bins, 1);
    end

    % Calculate significant wave height Hm0
    Hm0 = 4 * sqrt(m0);

end

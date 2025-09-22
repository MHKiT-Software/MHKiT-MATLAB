function m = frequency_moment(S, N, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the Nth frequency moment of the spectrum
%
% Parameters
% -----------
%     S: structure or numeric array
%         If structure:
%             S.spectrum   - Spectral density [m^2/Hz]
%             S.frequency  - Frequency [Hz]
%         If numeric:
%             S is spectral density array (vector or matrix)
%             varargin{1} must contain frequency vector
%
%     N: int
%         Moment (0 for 0th, 1 for 1st ....)
%
%     frequency_bins: vector (optional)
%         Bin widths for frequency of S. Required for unevenly sized bins
%
% Returns
% -------
%     m: double
%         Nth Frequency Moment indexed by S.columns
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    S
    N {mustBeNumeric, mustBeFinite, mustBeInteger}
end

arguments (Repeating)
    varargin
end

    % Extract spectrum and frequency
    if isstruct(S)
        spectrum = S.spectrum;
        frequency = S.frequency;
    elseif isnumeric(S)
        if nargin < 3
            error('When S is numeric, frequency vector must be provided as third argument');
        end
        spectrum = S;
        frequency = varargin{1};
        varargin(1) = [];
    else
        error('Input S must be a struct or numeric array');
    end

    % Standardize frequency, spectrum, and frequency bins
    if ~isempty(varargin)
        [frequency, spectrum, freq_bins] = standardize_wave_spectra_frequency(frequency, spectrum, varargin{1});
    else
        [frequency, spectrum, freq_bins] = standardize_wave_spectra_frequency(frequency, spectrum);
    end

    if isscalar(freq_bins)
        freq_bins = freq_bins(:);
    end

    % Calculate Nth moment: m_N = sum(f^N * S * df)
    m = sum((frequency.^N) .* spectrum .* freq_bins, 1);

end

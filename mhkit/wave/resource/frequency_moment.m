function m = frequency_moment(S, N, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the Nth frequency moment of the spectrum
%
% Computes: m_N = sum(f^N * S(f) * df) per Eq 8 in IEC 62600-101 Ed. 2.0 en 2024
%
% Parameters
% ------------
% S : struct or numeric
%   If struct:
%     S.spectrum : vector or matrix [m^2/Hz]
%       Spectral density
%     S.frequency : vector [Hz]
%       Frequency
%   If numeric: spectral density array, varargin{1} must be frequency
% N : integer
%   Moment order (0 for 0th, 1 for 1st, -1 for inverse, etc.)
% frequency_bins : vector (optional)
%   Bin widths for frequency of S. Required for unevenly sized bins.
%
% Returns
% ---------
% m : double or column vector
%   Nth frequency moment. Returns a column vector if S.spectrum is a matrix
%   (one value per column).
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
            error('MHKiT:frequency_moment:InvalidInput', ...
                'When S is numeric, frequency vector must be provided as third argument');
        end
        spectrum = S;
        frequency = varargin{1};
        varargin(1) = [];
    else
        error('MHKiT:frequency_moment:InvalidInput', 'Input S must be a struct or numeric array');
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
    m = m(:);
    mhkit_verify_column_vector_output(m, 'frequency_moment');

end

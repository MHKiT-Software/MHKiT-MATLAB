function Te = energy_period(S, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate wave energy period Te from wave spectra
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

arguments
    S
end

arguments (Repeating)
    varargin
end

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

    % Standardize frequency, spectrum, and frequency bins
    if ~isempty(varargin)
        [frequency, spectrum, freq_bins] = standardize_wave_spectra_frequency(frequency, spectrum, varargin{1});
    else
        [frequency, spectrum, freq_bins] = standardize_wave_spectra_frequency(frequency, spectrum);
    end

    % Calculate moments: m0 and m-1
    m0 = sum(spectrum .* freq_bins, 1);
    m_neg1 = sum((spectrum ./ frequency) .* freq_bins, 1);

    % Calculate energy period
    Te = m_neg1 ./ m0;

end

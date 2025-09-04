function Hm0 = significant_wave_height(S, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate significant wave height Hm0 from spectra
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
        error('Input S must be either a struct with fields .spectrum and .frequency, or a numeric array');
    end

    % Standardize frequency, spectrum, and frequency bins
    if ~isempty(varargin)
        [frequency, spectrum, freq_bins] = standardize_wave_spectra_frequency(frequency, spectrum, varargin{1});
    else
        [frequency, spectrum, freq_bins] = standardize_wave_spectra_frequency(frequency, spectrum);
    end

    % Calculate zeroth moment m0
    m0 = sum(spectrum .* freq_bins, 1);

    % Calculate significant wave height Hm0
    Hm0 = 4 * sqrt(m0);

end

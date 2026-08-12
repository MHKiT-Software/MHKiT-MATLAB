function [frequency, spectrum, freq_bins] = standardize_wave_spectra_frequency(frequency, spectrum, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standardize frequency and spectrum data for spectral calculations
%
% This utility function standardizes frequency vectors, spectra, and frequency bins
% following MHKiT-Python conventions for consistent numerical integration.
%
% Parameters
% ------------
%   frequency: vector
%       Frequency vector (Hz)
%   spectrum: vector or matrix  
%       Spectral density data
%   freq_bins: vector (optional)
%       Bin widths for frequency. If not provided, calculated from frequency differences
%
% Returns
% ---------
%   frequency: vector
%       Filtered frequency vector (> 1e-12 Hz), column format
%   spectrum: vector or matrix
%       Filtered spectrum data matching frequency
%   freq_bins: vector
%       Frequency bin widths, column format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    frequency {mustBeNumeric, mustBeFinite}
    spectrum {mustBeNumeric, mustBeFinite}
end

arguments (Repeating)
    varargin
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

    % Omit near-zero frequencies, following MHKiT-Python convention:
    % https://github.com/MHKiT-Software/MHKiT-Python/blob/6bad8fe4f2bd8a9bff66fb9607ed0900f09d0258/mhkit/wave/resource.py#L455-L457
    valid_idx = frequency > 1e-12;
    frequency = frequency(valid_idx);
    spectrum = spectrum(valid_idx, :);

    % Handle optional frequency bins input
    if ~isempty(varargin)
        freq_bins = varargin{1};
        if ~isscalar(freq_bins)
            freq_bins = freq_bins(:);
            if length(freq_bins) ~= length(valid_idx)
                error('Length of freq_bins must match original frequency vector');
            end
            freq_bins = freq_bins(valid_idx);
        else
            % Ensure scalar freq_bins is column vector for consistency
            freq_bins = freq_bins(:);
        end
    else
        % Bin widths computed after filtering, following MHKiT-Python convention:
        % https://github.com/MHKiT-Software/MHKiT-Python/blob/6bad8fe4f2bd8a9bff66fb9607ed0900f09d0258/mhkit/wave/resource.py#L460-L464
        % Computing them before filtering gives the wrong first bin width for non-uniform spacing.
        df = diff(frequency);
        freq_bins = [df(1); df(:)];  % Prepend first difference, ensure column vector
    end

end

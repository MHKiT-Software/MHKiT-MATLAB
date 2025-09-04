function [frequency, spectrum, freq_bins] = standardize_wave_spectra_frequency(frequency, spectrum, varargin)
%STANDARDIZE_WAVE_SPECTRA_FREQUENCY Standardize frequency and spectrum data for wave calculations
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

    % Handle optional frequency bins input
    if ~isempty(varargin)
        freq_bins = varargin{1};
    else
        % Calculate individual frequency bin widths to match MHKiT-Python implementation
        % MHKiT-Python uses: delta_f = f.diff() then prepends first difference
        % This creates a vector where each frequency has its own bin width,
        % which is critical for accurate numerical integration when frequencies
        % are not perfectly uniform (using mean(df) introduces systematic error)
        df = diff(frequency);
        freq_bins = [df(1); df(:)];  % Prepend first difference, ensure column vector
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

    % Filter out zero and near-zero frequencies to avoid division by zero
    % Following MHKiT-Python implementation: omit frequencies <= 1e-12
    valid_idx = frequency > 1e-12;
    frequency = frequency(valid_idx);
    spectrum = spectrum(valid_idx, :);
    
    % Handle frequency bins - filter if vector, keep scalar as-is
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

end
function Hm0 = significant_wave_height(S, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates significant wave height Hm0 from spectra
%
% Parameters
% ------------
%     S: structure or numeric array
%         If structure:
%             S.spectrum: Spectral density [m^2/Hz]
%             S.frequency: frequency [Hz]
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
%         Significant wave height [m] index by S.columns
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    S
end

arguments (Repeating)
    varargin
end

    % Calculate zeroth moment m0 using frequency_moment function
    if isstruct(S)
        % Pass struct directly
        m0 = frequency_moment(S, 0, varargin{:});
    elseif isnumeric(S)
        % Pass numeric spectrum with frequency vector
        if nargin < 2
            error('When S is numeric, frequency vector must be provided as second argument');
        end
        m0 = frequency_moment(S, 0, varargin{:});
    else
        error('Input S must be either a struct with fields .spectrum and .frequency, or a numeric array');
    end

    % Calculate significant wave height Hm0
    Hm0 = 4 * sqrt(m0);

end

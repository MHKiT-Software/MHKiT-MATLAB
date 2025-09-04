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

    % Calculate moments using frequency_moment function
    if isstruct(S)
        % Pass struct directly
        m0 = frequency_moment(S, 0, varargin{:});
        m_neg1 = frequency_moment(S, -1, varargin{:});
    elseif isnumeric(S)
        % Pass numeric spectrum with frequency vector
        if nargin < 2
            error('When S is numeric, frequency vector must be provided as second argument');
        end
        m0 = frequency_moment(S, 0, varargin{:});
        m_neg1 = frequency_moment(S, -1, varargin{:});
    else
        error('Input S must be a struct or numeric array');
    end

    % Calculate energy period
    Te = m_neg1 ./ m0;

end

function Tp = peak_period(S)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates wave peak period from spectra
%
% Peak period is the inverse of the frequency at which the spectrum
% has maximum energy (Eq 14 in IEC 62600-101 Ed. 2.0 en 2024).
%
% Parameters
% ------------
% S : struct
%   Wave spectrum structure:
%     S.spectrum : vector or matrix [m^2/Hz]
%       Spectral density
%     S.frequency : vector [Hz]
%       Frequency
%     S.type : string (optional)
%       Spectra type description
%
% Returns
% ---------
% Tp : double or vector [s]
%   Wave peak period. Returns a vector if S.spectrum is a matrix
%   (one value per column).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    S struct
end

arguments (Output)
    Tp {mustBeNumeric}
end

% Validate input structure
if ~isfield(S, 'spectrum') || ~isfield(S, 'frequency')
    error('MHKiT:peak_period:InvalidInput', ...
        'S must be a structure with spectrum and frequency fields');
end

spectrum = S.spectrum;
frequency = S.frequency(:);

% Handle both vector and matrix spectra
if isvector(spectrum)
    spectrum = spectrum(:);
    [~, idx] = max(spectrum);
    fp = frequency(idx);
else
    % Matrix case: find max along first dimension (frequency)
    [~, idx] = max(spectrum, [], 1);
    fp = frequency(idx);
end

% Eq 14 in IEC 62600-101 Ed. 2.0 en 2024
Tp = 1 ./ fp;

end

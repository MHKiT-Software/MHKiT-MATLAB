function fmax_out = fmax_warning(fn, fmax)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check maximum frequency limit isn't greater than the Nyquist frequency
%
% Parameters
% ------------
%   fn: double (scalar)
%       The Nyquist frequency [Hz]
%   fmax: double (scalar)
%       The maximum frequency limit [Hz]
%
% Returns
% ---------
%   fmax_out: double (scalar)
%       The adjusted maximum frequency limit, ensuring it does not exceed
%       the Nyquist frequency [Hz]
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    fn {mustBeNumeric}
    fmax {mustBeNumeric}
end

arguments (Output)
    fmax_out {mustBeNumeric}
end

% Validate inputs are positive
if fn <= 0
    error('MHKiT:acoustics:fmax_warning:InvalidNyquist', ...
        'Nyquist frequency must be positive. Got fn = %.3f', fn);
end

if fmax <= 0
    error('MHKiT:acoustics:fmax_warning:InvalidFmax', ...
        'Maximum frequency must be positive. Got fmax = %.3f', fmax);
end

if fmax > fn
    warning('MHKiT:acoustics:fmax_warning:FrequencyExceedsNyquist', ...
        'fmax = %.1f is greater than Nyquist frequency. Setting fmax = %.1f', fmax, fn);
    fmax = fn;
end

fmax_out = fmax;

end
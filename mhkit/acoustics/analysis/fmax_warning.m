function fmax_out = fmax_warning(fn, fmax)
% 
% Checks that the maximum frequency limit isn't greater than the Nyquist frequency.
% 
% Parameters
% ----------
% fn: double
%     The Nyquist frequency in Hz.
% fmax: double
%     The maximum frequency limit in Hz.
% 
% Returns
% -------
% fmax: double
%     The adjusted maximum frequency limit, ensuring it does not exceed 
%     the Nyquist frequency.


arguments (Input)
    fn {mustBeNumeric}
    fmax {mustBeNumeric}
end
    
arguments (Output)
    fmax_out {mustBeNumeric}
end

if fmax > fn
    warning('fmax = %.1f is greater than Nyquist frequency. Setting fmax = %.1f', fmax, fn)
    fmax = fn;
end

fmax_out = fmax;

end
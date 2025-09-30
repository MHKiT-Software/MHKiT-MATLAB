function win = mhkit_get_window(window_type, n_fft)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Dispatch window function for signal processing applications
%
% Parameters
% ------------
%   window_type: char or double
%       Window type {'hann', 'hamming', 'rectwin'} or custom window vector
%   n_fft: double
%       Length of the window (number of samples)
%
% Returns
% ---------
%   win: double
%       Window values as column vector [n_fft x 1]
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    window_type
    n_fft (1,1) {mustBeNumeric, mustBePositive, mustBeInteger}
end

% Create window based on type
if ischar(window_type) || isstring(window_type)
    switch lower(char(window_type))
        case 'hann'
            win = mhkit_window_hann(n_fft);
        case 'hamming'
            win = mhkit_window_hamming(n_fft);
        case 'rectwin'
            win = mhkit_window_rectwin(n_fft);
        otherwise
            error('MHKiT:mhkit_get_window:InvalidWindow', 'Invalid window type: %s', char(window_type));
    end
elseif isnumeric(window_type)
    win = double(window_type(:));
    if length(win) ~= n_fft
        error('MHKiT:mhkit_get_window:InvalidWindow', 'Custom window length must equal n_fft');
    end
else
    error('MHKiT:mhkit_get_window:InvalidWindow', 'Window type must be string or numeric vector');
end

end

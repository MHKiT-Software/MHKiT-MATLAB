function win = mhkit_window_hann(n_fft)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Generate Hann window for signal processing applications
%
% Parameters
% ------------
%   n_fft: double
%       Length of the window (number of samples)
%
% Returns
% ---------
%   win: double
%       Hann window values as column vector [n_fft x 1]
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    n_fft (1,1) {mustBeNumeric, mustBePositive}
end

n_fft = double(n_fft);

% Create sample indices 0 to N-1)
n = (0:n_fft-1)';

% Compute Hann window
win = 0.5 - 0.5 * cos(2*pi*n/(n_fft-1));

end

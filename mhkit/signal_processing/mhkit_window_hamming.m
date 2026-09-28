function win = mhkit_window_hamming(n_fft)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Generate Hamming window for signal processing applications
%
% Parameters
% ------------
%   n_fft: double
%       Length of the window (number of samples)
%
% Returns
% ---------
%   win: double
%       Hamming window values as column vector [n_fft x 1]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    n_fft (1,1) {mustBeNumeric, mustBePositive, mustBeInteger}
end

n_fft = double(n_fft);

% Create sample indices (0 to N-1)
n = (0:n_fft-1)';

% Compute Hamming window
win = 0.54 - 0.46 * cos(2*pi*n/(n_fft-1));

end

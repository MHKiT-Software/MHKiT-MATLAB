function win = mhkit_window_rectwin(n_fft)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Generate rectangular window for signal processing applications
%
% Parameters
% ------------
%   n_fft: double
%       Length of the window (number of samples)
%
% Returns
% ---------
%   win: double
%       Rectangular window values as column vector [n_fft x 1]
%       All values are 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    n_fft (1,1) {mustBeNumeric, mustBePositive, mustBeInteger}
end

n_fft = double(n_fft);

% Create rectangular window (all ones)
win = ones(n_fft, 1);

end

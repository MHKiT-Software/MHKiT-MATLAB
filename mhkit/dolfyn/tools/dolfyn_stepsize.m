function [step, n_segments, n_fft_used] = dolfyn_stepsize(l, n_fft, n_segments, step)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Calculate step size for overlapping FFT segments
%
%     Use same algorithm as MHKiT-Python _stepsize function
%
% Parameters
% ------------
%   l : double
%       Length of input data array [samples]
%   n_fft : double  
%       Desired FFT window size [samples]
%   n_segments : double, optional
%       Number of overlapping segments (default: auto-calculated)
%   step : double, optional
%       Step size between segments (default: auto-calculated) [samples]
%
% Returns
% ---------
%   step : double
%       Step size for segment overlap [samples]
%   n_segments : double
%       Number of overlapping segments [dimensionless]
%   n_fft_used : double
%       Actual FFT size used (may be reduced if l < n_fft) [samples]
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    l (1,1) {mustBeNumeric, mustBePositive}
    n_fft (1,1) {mustBeNumeric, mustBePositive}
    n_segments {mustBeNumeric, mustBeNonnegative} = []
    step {mustBeNumeric, mustBeNonnegative} = []
end

% Additional input validation for non-empty values
if ~isempty(n_segments) && n_segments < 1
    error('MHKiT:dolfyn_stepsize:InvalidSegments', ...
          'Number of segments must be positive');
end

% Validate that FFT size is not larger than data length
if n_fft > l
    warning('MHKiT:dolfyn_stepsize:FFTSizeReduced', ...
            'FFT size (%d) larger than data length (%d), reducing to data length', ...
            n_fft, l);
end

% Ensure inputs are integers
l = double(l);
n_fft = double(n_fft);

% MHKiT-Python: if l < nfft, adjust nfft (handled in calling function)
n_fft_used = n_fft;
if l < n_fft
    n_fft_used = l;
end

% MHKiT-Python _stepsize algorithm
if isempty(n_segments) && isempty(step)
    % Auto-calculate both (most common case)
    if l == n_fft_used
        % Special case: single segment
        step = 0;
        n_segments = 1;
    else
        % MHKiT-Python: nens = int(2.0 * l / nfft)
        n_segments = floor(2.0 * l / n_fft_used);
        % MHKiT-Python: step = int((l - nfft) / (nens - 1))
        step = floor((l - n_fft_used) / (n_segments - 1));
    end
elseif isempty(step)
    % n_segments specified, calculate step
    % MHKiT-Python: step = int((l - nfft) / (nens - 1))
    step = floor((l - n_fft_used) / (n_segments - 1));
else
    % step specified, calculate n_segments
    % MHKiT-Python: nens = int((l - nfft) / step + 1)
    n_segments = floor((l - n_fft_used) / step + 1);
end

% Ensure outputs are integers
step = double(step);
n_segments = double(n_segments);
n_fft_used = double(n_fft_used);

end

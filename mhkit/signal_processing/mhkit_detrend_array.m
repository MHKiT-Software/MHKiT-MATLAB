function detrended = mhkit_detrend_array(data)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Detrend array by removing mean and linear trend
%
% Parameters
% ------------
%   data: double
%       Input data vector to detrend
%
% Returns
% ---------
%   detrended: double
%       Detrended data with same dimensions as input
%       Mean and linear trend removed
%
% Key Equations
% -------------
% 1. Mean removal:
%    data_centered = data - mean(data)
%
% 2. Linear trend removal:
%    slope = mean(x * data_centered) / mean(x^2)
%    detrended = data_centered - slope * x_centered
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    data {mustBeNumeric}
end

% Validate input is not empty
if isempty(data)
    error('MHKiT:mhkit_detrend_array:EmptyInput', 'Input data cannot be empty');
end

data = double(data(:));

% Handle edge cases that cause matrix dimension errors
if length(data) < 2
    detrended = data - mhkit_nanmean(data);  % Just remove mean
    return;
end

% Check for all NaN or constant data
valid_count = sum(~isnan(data));
if valid_count < 2
    detrended = data - mhkit_nanmean(data);  % Just remove mean
    return;
end

% Step 1: Remove mean first (MHKiT-Python: arr -= np.nanmean(arr))
data_centered = data - mhkit_nanmean(data);

% Step 2: Create and center index vector (MHKiT-Python: x -= np.nanmean(x))
x = (0:length(data)-1)';
x_centered = x - mhkit_nanmean(x);

% Step 3: Compute slope (MHKiT-Python: b = np.nanmean(x*arr) / np.nanmean(x**2))
x_sq_mean = mhkit_nanmean(x_centered.^2);
if abs(x_sq_mean) > eps  % Avoid division by zero
    xy_mean = mhkit_nanmean(x_centered .* data_centered);
    slope = xy_mean / x_sq_mean;
    % Step 4: Remove linear trend
    detrended = data_centered - slope * x_centered;
else
    % No trend to remove (constant or single point)
    detrended = data_centered;
end

end

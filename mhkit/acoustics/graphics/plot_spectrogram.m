function plot_spectrogram(spsdl, options)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot spectrogram of the sound pressure spectral density level
%
% Parameters
% ------------
%   spsdl: struct
%       spsdl.data : Spectral density level data [dB rel 1 uPa^2/Hz]
%       spsdl.freq : Frequency vector [Hz]
%       spsdl.time : Time vector
%       spsdl.units : Data units
%   options: name-value pairs (optional)
%       fmin: float - Lower frequency limit [Hz] (default: 10)
%       fmax: float - Upper frequency limit [Hz] (default: 100000)
%       cm: string - Colormap name (default: 'hot')
%       cmin: float - Minimum colorbar value (default: [])
%       cmax: float - Maximum colorbar value (default: [])
%       use_smoothing: logical - Apply median smoothing (default: true)
%       smoothing_stride: integer - Smoothing kernel size (default: 3)
%       plot_method: string - 'pcolor' (slowest), 'surf' (medium), or 'imagesc' (fastest) (default: 'imagesc')
%
% Returns
% ---------
%   This function creates a figure but does not return a value
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    spsdl struct
    options.fmin {mustBeNumeric} = 10
    options.fmax {mustBeNumeric} = 100000
    options.cm {mustBeText} = 'hot'
    options.cmin {mustBeNumeric} = []
    options.cmax {mustBeNumeric} = []
    options.use_smoothing {mustBeNumericOrLogical} = true
    options.smoothing_stride {mustBeInteger, mustBePositive} = 3
    options.plot_method {mustBeText} = 'imagesc'
end

% Validate spsdl structure
validate_spsdl_struct(spsdl, 'plot_spectrogram', ...
    'required_fields', {{'data', 'freq', 'time', 'units'}}, ...
    'check_dimensions', false);

% Validate plot method
if ~ismember(options.plot_method, {'pcolor', 'surf', 'imagesc'})
    error('MHKiT:acoustics:plot_spectrogram:InvalidInput', ...
        'plot_method must be pcolor, surf, or imagesc');
end

% Check fmax
fn = max(spsdl.freq);
fmax = fmax_warning(fn, options.fmax);

% Select frequency range
idxf = find(spsdl.freq>=options.fmin & spsdl.freq<=fmax);
data_subset = spsdl.data(idxf,:);

% Apply smoothing if requested
if options.use_smoothing
    data_subset = movmedian(movmedian(data_subset, options.smoothing_stride, 1), options.smoothing_stride, 2);
end

figure;

% Plot using specified method
switch options.plot_method
    case 'imagesc'
        imagesc(spsdl.time, spsdl.freq(idxf), data_subset);
        set(gca, 'YDir', 'normal');
    case 'surf'
        surf(spsdl.time, spsdl.freq(idxf), data_subset, 'EdgeColor', 'none');
        view(0, 90);
    case 'pcolor'
        pcolor(spsdl.time, spsdl.freq(idxf), data_subset);
        shading interp;
end

ylim([options.fmin fmax]);
set(gca, 'YScale', 'log');
colormap(options.cm);
cb = colorbar();

% Set colorbar limits based on what's defined
if ~isempty(options.cmin) && ~isempty(options.cmax)
    % Both limits defined
    clim([options.cmin options.cmax]);
elseif ~isempty(options.cmin)
    % Only minimum defined
    current_limits = clim;
    clim([options.cmin current_limits(2)]);
elseif ~isempty(options.cmax)
    % Only maximum defined
    current_limits = clim;
    clim([current_limits(1) options.cmax]);
end
% If neither defined, use MATLAB default behavior

xlabel('Time');
ylabel('Frequency [Hz]');
cb.Label.String = spsdl.units;

end


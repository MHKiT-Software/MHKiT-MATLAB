function ds_out = calculate_doppler_noise_level(ds, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate bias due to Doppler noise using the noise floor of velocity spectra.
%
% This function estimates the Doppler noise level from the high-frequency
% region of power spectral density data. The noise level is calculated by
% examining the spectral density at frequencies approaching the Nyquist limit.
%
% Parameters
% ------------
%   ds: structure
%       ADCP dataset structure containing power spectral density data
%   psd_field: string
%       Name of the PSD field in the dataset
%       Default: 'psd'
%   pct_fN: numeric
%       Percent of Nyquist frequency to define characteristic frequency [0-1]
%       Default: 0.8 (80% of Nyquist frequency)
%   field_name: string
%       Name of output field in dataset structure
%       Default: 'noise_level'
%
% Returns
% ---------
%   ds_out: structure
%       Input dataset with added Doppler noise level field:
%           ds_out.noise_level.data : Noise level values [m/s]
%           ds_out.noise_level.dims : dimension names
%           ds_out.noise_level.coords : coordinate information
%           ds_out.noise_level.units : "m s-1"
%           ds_out.noise_level.long_name : "Doppler Noise Level"
%
% Example
% -------
% % Calculate noise from existing PSD
% ds_noise = calculate_doppler_noise_level(ds_with_psd, 'pct_fN', 0.8);
%
% References
% ----------
% Richard, Jean-Baptiste, et al. "Method for identification of Doppler noise
% levels in turbulent flow measurements dedicated to tidal energy." International
% Journal of Marine Energy 3 (2013): 52-64.
%
% Thiébaut, Maxime, et al. "Investigating the flow dynamics and turbulence at a
% tidal-stream energy site in a highly energetic estuary." Renewable Energy 195
% (2022): 252-262.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    arguments
        ds
        options.psd_field = 'psd'
        options.pct_fN = 0.8
        options.field_name = 'noise_level'
    end
    
    % Validate input dataset structure
    if ~isstruct(ds)
        error('MHKiT:dolfyn:calculate_doppler_noise_level: Input ds must be a structure');
    end
    
    % Check for required PSD field
    if ~isfield(ds, options.psd_field)
        error('MHKiT:dolfyn:calculate_doppler_noise_level: Dataset must contain PSD field: %s', options.psd_field);
    end
    
    psd_data = ds.(options.psd_field);
    
    % Validate PSD structure
    if ~isfield(psd_data, 'data') || ~isfield(psd_data, 'dims') || ~isfield(psd_data, 'coords')
        error('MHKiT:dolfyn:calculate_doppler_noise_level: PSD field must contain data, dims, and coords');
    end
    
    % Validate pct_fN parameter
    if options.pct_fN < 0 || options.pct_fN > 1
        error('MHKiT:dolfyn:calculate_doppler_noise_level: pct_fN must be between 0 and 1');
    end
    
    % Get frequency coordinate
    if ~isfield(psd_data.coords, 'freq')
        error('MHKiT:dolfyn:calculate_doppler_noise_level: PSD data must contain freq coordinate');
    end
    
    freq_data = psd_data.coords.freq;
    psd_values = psd_data.data;
    
    % Determine frequency units and calculate characteristic frequency
    % Characteristic frequency set to 80% of Nyquist frequency
    if isfield(psd_data.coords, 'freq') && isfield(ds, 'attrs') && isfield(ds.attrs, 'fs')
        fs = ds.attrs.fs;
        fN = fs / 2;  % Nyquist frequency
    else
        % Estimate Nyquist from frequency data
        fN = max(freq_data);
    end
    
    fc = options.pct_fN * fN;
    
    % Determine frequency units for proper scaling
    if isfield(psd_data, 'attrs') && isfield(psd_data.attrs, 'units') && contains(lower(string(psd_data.attrs.units)), 'hz')
        freq_range_mask = (freq_data >= fc) & (freq_data <= fN);
    else
        % Assume rad/s units
        fc_rad = 2 * pi * fc;
        fN_rad = 2 * pi * fN;
        freq_range_mask = (freq_data >= fc_rad) & (freq_data <= fN_rad);
    end
    
    % Check if we have valid frequency range
    if sum(freq_range_mask) == 0
        error('MHKiT:dolfyn:calculate_doppler_noise_level: No frequency points found in range [%.3f, %.3f]', fc, fN);
    end
    
    % Find frequency dimension in PSD data
    freq_dim = find(strcmp(psd_data.dims, 'freq'));
    if isempty(freq_dim)
        error('MHKiT:dolfyn:calculate_doppler_noise_level: PSD data must have freq dimension');
    end
    
    % Calculate noise floor in high-frequency region
    % Extract PSD values in the noise floor region
    if freq_dim == 1
        psd_noise_region = psd_values(freq_range_mask, :);
        freq_noise_region = freq_data(freq_range_mask);
    elseif freq_dim == 2
        psd_noise_region = psd_values(:, freq_range_mask);
        freq_noise_region = freq_data(freq_range_mask);
    elseif freq_dim == 3
        psd_noise_region = psd_values(:, :, freq_range_mask);
        freq_noise_region = freq_data(freq_range_mask);
    else
        error('MHKiT:dolfyn:calculate_doppler_noise_level: Frequency dimension must be 1, 2, or 3');
    end
    
    % Calculate N2 = PSD * freq for noise estimation
    if freq_dim == 1
        N2 = psd_noise_region .* freq_noise_region;
    elseif freq_dim == 2
        N2 = psd_noise_region .* freq_noise_region';
    else  % freq_dim == 3
        % For 3D case [range, dir, freq], broadcast frequency along 3rd dimension
        N2 = psd_noise_region .* reshape(freq_noise_region, [1, 1, length(freq_noise_region)]);
    end
    
    % Calculate noise level: sqrt(mean(N2))
    noise_level = sqrt(mean(N2, freq_dim, 'omitnan'));
    
    % Create output coordinates (remove freq dimension)
    output_coords = psd_data.coords;
    output_coords = rmfield(output_coords, 'freq');
    
    % Create output dimensions (remove freq dimension)
    output_dims = psd_data.dims(~strcmp(psd_data.dims, 'freq'));
    
    % Create output structure
    ds_out = ds;
    ds_out.(options.field_name) = struct();
    ds_out.(options.field_name).data = single(noise_level);
    ds_out.(options.field_name).dims = output_dims;
    ds_out.(options.field_name).coords = output_coords;
    ds_out.(options.field_name).units = "m s-1";
    ds_out.(options.field_name).long_name = "Doppler Noise Level";
    ds_out.(options.field_name).description = sprintf(...
        'Doppler noise level calculated from PSD white noise at %.1f%% of Nyquist frequency', ...
        options.pct_fN * 100);

end

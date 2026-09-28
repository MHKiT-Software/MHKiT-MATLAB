function ds_out = calculate_dissipation_rate_LT83(ds, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate turbulent kinetic energy dissipation rate from PSD using Lumley-Terray 1983 method.
%
% This function calculates dissipation rate from power spectral density data
% using the inertial subrange method. The calculation applies the Kolmogorov
% -5/3 law to relate spectral density to dissipation rate.
%
% Parameters
% ------------
%   ds: structure
%       ADCP dataset structure containing PSD and velocity magnitude data
%   psd_field: string
%       Name of the PSD field in the dataset
%       Default: 'psd'
%   U_mag_field: string
%       Name of the velocity magnitude field in the dataset
%       Default: 'U_mag'
%   freq_range: array [1x2]
%       Frequency range for inertial subrange [f_min, f_max] in Hz or rad/s
%       Default: [0.2, 0.4]
%   noise: numeric or structure
%       Doppler noise level in same units as velocity [m/s]
%       Can be scalar, array, or dataset field structure
%       Default: 0 (no noise correction)
%   field_name: string
%       Name of output field in dataset structure
%       Default: 'dissipation_rate'
%
% Returns
% ---------
%   ds_out: structure
%       Input dataset with added dissipation rate field:
%           ds_out.dissipation_rate.data : Dissipation rate values [m² s⁻³]
%           ds_out.dissipation_rate.dims : dimension names
%           ds_out.dissipation_rate.coords : coordinate information
%           ds_out.dissipation_rate.units : "m2 s-3"
%           ds_out.dissipation_rate.long_name : "TKE Dissipation Rate"
%
% References
% ----------
% LT83 : Lumley and Terray, "Kinematics of turbulence convected
% by a random wave field". JPO, 1983, vol13, pp2000-2007.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    arguments
        ds
        options.psd_field = 'psd'
        options.U_mag_field = 'U_mag'
        options.freq_range = [0.2, 0.4]
        options.noise = 0
        options.field_name = 'dissipation_rate'
    end
    
    % Validate input dataset structure
    if ~isstruct(ds)
        error('mhkit:dolfyn:calculate_dissipation_rate_LT83: Input ds must be a structure');
    end
    
    % Check for required PSD field
    if ~isfield(ds, options.psd_field)
        error('mhkit:dolfyn:calculate_dissipation_rate_LT83: Dataset must contain PSD field: %s', options.psd_field);
    end
    
    % Check for required velocity magnitude field
    if ~isfield(ds, options.U_mag_field)
        error('mhkit:dolfyn:calculate_dissipation_rate_LT83: Dataset must contain velocity magnitude field: %s', options.U_mag_field);
    end
    
    psd_data = ds.(options.psd_field);
    U_mag_data = ds.(options.U_mag_field);
    
    % Validate PSD structure
    if ~isfield(psd_data, 'data') || ~isfield(psd_data, 'dims') || ~isfield(psd_data, 'coords')
        error('mhkit:dolfyn:calculate_dissipation_rate_LT83: PSD field must contain data, dims, and coords');
    end
    
    % Validate U_mag structure
    if ~isfield(U_mag_data, 'data')
        error('mhkit:dolfyn:calculate_dissipation_rate_LT83: U_mag field must contain data');
    end
    
    % Validate frequency range
    if length(options.freq_range) ~= 2 || options.freq_range(1) >= options.freq_range(2)
        error('mhkit:dolfyn:calculate_dissipation_rate_LT83: freq_range must be [f_min, f_max] with f_min < f_max');
    end
    
    % Get frequency coordinate
    if ~isfield(psd_data.coords, 'freq')
        error('mhkit:dolfyn:calculate_dissipation_rate_LT83: PSD data must contain freq coordinate');
    end
    
    freq_data = psd_data.coords.freq;
    psd_values = psd_data.data;
    U_mag_values = U_mag_data.data;
    
    % Handle both 2D and 3D PSD data
    % 2D: [time, freq] - legacy format
    % 3D: [range, dir, freq] - current ADCP format
    psd_ndims = ndims(psd_values);
    
    if psd_ndims == 2
        process_3d = false;
    elseif psd_ndims == 3
        process_3d = true;
        [n_range, n_dir, n_freq] = size(psd_values);
    else
        error('mhkit:dolfyn:calculate_dissipation_rate_LT83: PSD data should be 2D [time, freq] or 3D [range, dir, freq]');
    end
    
    % Validate and handle U_mag dimensions
    if process_3d
        % For 3D PSD, U_mag should be [time, range] or [range, time]
        if ndims(U_mag_values) ~= 2
            error('mhkit:dolfyn:calculate_dissipation_rate_LT83: For 3D PSD, U_mag data should be 2-dimensional (time, range)');
        end
        
        % Check if dimensions match - U_mag could be [time, range] or [range, time]
        [u_dim1, u_dim2] = size(U_mag_values);
        if u_dim2 == n_range && u_dim1 ~= n_range
            % U_mag is [time, range]
            U_mag_values = U_mag_values';  % Transpose to [range, time]
        else
            error('mhkit:dolfyn:calculate_dissipation_rate_LT83: U_mag dimensions [%d, %d] do not match PSD range dimension %d', ...
                u_dim1, u_dim2, n_range);
        end
    else
        % Legacy 2D case
        if ndims(U_mag_values) ~= 2 || min(size(U_mag_values)) ~= 1
            error('mhkit:dolfyn:calculate_dissipation_rate_LT83: U_mag data should be 1-dimensional (time)');
        end
        
        % Ensure U_mag is column vector to match PSD time dimension
        U_mag_values = U_mag_values(:);
        
        % Check dimension compatibility
        if size(psd_values, 1) ~= length(U_mag_values)
            error('mhkit:dolfyn:calculate_dissipation_rate_LT83: PSD and U_mag must have same time dimension');
        end
    end
    
    % Handle noise parameter for both 2D and 3D cases
    if isstruct(options.noise) && isfield(options.noise, 'data')
        noise_values = options.noise.data;
    else
        noise_values = options.noise;
    end
    
    % Validate and reshape noise dimensions
    if process_3d
        % For 3D case, noise should match [range, dir] or be expandable
        if isscalar(noise_values)
            noise_values = noise_values * ones(n_range, n_dir);
        else
            noise_size = size(noise_values);
            if isequal(noise_size, [n_range, n_dir])
                % Already correct size
            elseif length(noise_values) == n_range
                % Expand to [range, dir]
                noise_values = repmat(noise_values(:), 1, n_dir);
            else
                error('mhkit:dolfyn:calculate_dissipation_rate_LT83: For 3D PSD, noise must be scalar, [range] or [range, dir]');
            end
        end
    else
        % Legacy 2D case
        if isscalar(noise_values)
            noise_values = noise_values * ones(size(U_mag_values));
        else
            noise_values = noise_values(:);
        end
        
        if length(noise_values) ~= length(U_mag_values)
            error('mhkit:dolfyn:calculate_dissipation_rate_LT83: Noise must have same length as U_mag time dimension');
        end
    end
    
    % Get sampling frequency for noise normalization
    if isfield(ds, 'attrs') && isfield(ds.attrs, 'fs')
        fs = ds.attrs.fs;
    else
        error('mhkit:dolfyn:calculate_dissipation_rate_LT83: Sampling frequency fs not found in ds.attrs');
    end
    
    % Select frequency range for inertial subrange
    freq_mask = (freq_data >= options.freq_range(1)) & (freq_data <= options.freq_range(2));
    
    if sum(freq_mask) == 0
        error('mhkit:dolfyn:calculate_dissipation_rate_LT83: No frequency points found in range [%.3f, %.3f]', ...
            options.freq_range(1), options.freq_range(2));
    end

    % Extract frequency range
    freq_inertial = freq_data(freq_mask);
    
    % Check frequency units and convert velocity accordingly
    if isfield(psd_data, 'units') && contains(lower(string(psd_data.units)), 'hz')
        % Frequency in Hz, convert velocity for rad/s equivalent
        velocity_scaling = 1 / (2 * pi);
    else
        % Assume rad/s units
        velocity_scaling = 1;
    end
    
    % Kolmogorov constant for single velocity component
    alpha = 0.5;
    
    % Calculate dissipation rate using Lumley-Terray 1983 method
    % ε = [(S(f) * f^(5/3) / α)^(3/2)] / U
    
    if process_3d
        % Process 3D PSD data [range, dir, freq]
        dissipation_rate = zeros(n_range, n_dir, 'single');
        processed_count = 0;
        
        for r = 1:n_range
            for d = 1:n_dir
                % Extract PSD for this range bin and component
                psd_series = squeeze(psd_values(r, d, :));  % [freq]
                
                % Skip if insufficient valid data
                valid_psd = psd_series(freq_mask);
                if sum(~isnan(valid_psd) & valid_psd > 0) < length(freq_inertial)/2
                    dissipation_rate(r, d) = NaN;
                    continue;
                end
                
                % Apply noise correction
                if any(noise_values ~= 0)
                    noise_correction = (noise_values(r, d)^2) / (fs / 2);
                    valid_psd = valid_psd - noise_correction;

                    % Handle negative/zero values after noise correction
                    negative_mask = valid_psd <= 0;
                    if any(negative_mask)
                        positive_values = valid_psd(valid_psd > 0);
                        if ~isempty(positive_values)
                            replacement_value = min(abs(positive_values)) / 100;
                        else
                            replacement_value = 1e-10;  % Fallback value
                        end
                        valid_psd(negative_mask) = replacement_value;
                    end
                end
                
                % Get corresponding velocity magnitude
                if size(U_mag_values, 1) >= r
                    U_effective = mean(U_mag_values(r, :), 'omitnan') * velocity_scaling;
                else
                    U_effective = NaN;
                end
                
                if U_effective <= 0 || isnan(U_effective)
                    dissipation_rate(r, d) = NaN;
                    continue;
                end
                
                % Compute S(f) * f^(5/3)
                freq_factor = freq_inertial.^(5/3);
                spectral_integral = valid_psd .* freq_factor;
                
                % Average over frequency range
                mean_spectral = mean(spectral_integral, 'omitnan');
                
                % Apply power law and normalize by velocity
                eps_value = (mean_spectral / alpha).^(3/2) / U_effective;
                
                dissipation_rate(r, d) = max(0, eps_value);  % Ensure positive
                processed_count = processed_count + 1;
            end
        end

    else
        % Legacy 2D processing [time, freq]
        psd_corrected = psd_values;
        
        % Apply noise correction to PSD
        if any(noise_values ~= 0)
            noise_correction = (noise_values.^2) / (fs / 2);
            psd_corrected = psd_corrected - noise_correction;
            
            % Set minimum threshold to prevent negative/zero PSDs
            min_psd = min(abs(psd_corrected(:))) / 100;
            psd_corrected(psd_corrected <= 0) = min_psd;
        end
        
        U_effective = U_mag_values * velocity_scaling;
        
        % Extract PSD values in frequency range  
        psd_inertial = psd_corrected(:, freq_mask);
        
        % Compute S(f) * f^(5/3) for each time step
        freq_factor = freq_inertial.^(5/3);
        spectral_integral = psd_inertial .* freq_factor';
        
        % Average over frequency range
        mean_spectral = mean(spectral_integral, 2, 'omitnan');
        
        % Apply power law and normalize by velocity
        dissipation_rate = (mean_spectral / alpha).^(3/2) ./ U_effective;
        
        % Handle any invalid results
        dissipation_rate(U_effective <= 0) = NaN;
        dissipation_rate(dissipation_rate < 0) = NaN;
    end
    
    % Create output structure based on data format
    if process_3d
        % For 3D case, create coordinates and dimensions for [range, dir]
        output_coords = struct();

        % Get range coordinates from U_mag field (which should have both time and range)
        if isfield(U_mag_data, 'coords') && isfield(U_mag_data.coords, 'range')
            output_coords.range = U_mag_data.coords.range;
        else
            % Fallback: create range indices
            output_coords.range = 1:n_range;
        end

        % For direction coordinate, use default labels since PSD doesn't have this info
        output_coords.dir = 1:n_dir;
        output_dims = {'range', 'dir'};
        
        % Calculate summary statistics
        valid_estimates = sum(~isnan(dissipation_rate(:)));
        total_estimates = numel(dissipation_rate);
    else
        % Legacy 2D case - use U_mag structure
        output_coords = U_mag_data.coords;
        output_dims = U_mag_data.dims;
    end
    
    % Create output structure
    ds_out = ds;
    ds_out.(options.field_name) = struct();
    ds_out.(options.field_name).data = single(dissipation_rate);
    ds_out.(options.field_name).dims = output_dims;
    ds_out.(options.field_name).coords = output_coords;
    ds_out.(options.field_name).units = "m2 s-3";
    ds_out.(options.field_name).long_name = "TKE Dissipation Rate";
    ds_out.(options.field_name).description = sprintf(...
        'TKE dissipation rate calculated using Lumley-Terray 1983 method in freq range [%.3f, %.3f]', ...
        options.freq_range(1), options.freq_range(2));
        
    % Add method metadata
    ds_out.(options.field_name).method = 'Lumley-Terray 1983';
    ds_out.(options.field_name).kolmogorov_constant = alpha;
    ds_out.(options.field_name).frequency_range = options.freq_range;
    ds_out.(options.field_name).noise_corrected = any(noise_values(:) ~= 0);
end

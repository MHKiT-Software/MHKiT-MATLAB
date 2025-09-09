function ds_out = calculate_turbulence_intensity(ds, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate noise-corrected turbulence intensity from ADCP velocity data.
% 
% This function calculates turbulence intensity from ensemble-averaged ADCP data.
% Input data should be pre-averaged using average_by_dimension() function.
% TI is calculated as the ratio of velocity fluctuation RMS to mean velocity.
%
% Parameters
% ------------
%   ds: structure
%       ADCP dataset structure containing ensemble-averaged velocity data
%   noise: numeric or array
%       Doppler noise level in same units as velocity [m/s]
%       Default: 0 (no noise correction)
%   thresh: numeric  
%       Velocity threshold below which TI will not be calculated [m/s]
%       Default: 0 (no threshold)
%   field_name: string
%       Name of output field in dataset structure
%       Default: 'turbulence_intensity'
%
% Returns
% ---------
%   ds_out: structure
%       Input dataset with added turbulence intensity field:
%           ds_out.turbulence_intensity.data : TI values [dimensionless]
%           ds_out.turbulence_intensity.dims : dimension names  
%           ds_out.turbulence_intensity.coords : coordinate information
%           ds_out.turbulence_intensity.units : "1"
%           ds_out.turbulence_intensity.long_name : "Turbulence Intensity"
%
% Example
% -------
% % First, ensemble average the raw data
% ds_avg = average_by_dimension(ds_raw, 300, 'time');  % 300-sample ensembles
% % Then calculate turbulence intensity  
% ds_avg = calculate_turbulence_intensity(ds_avg, 'noise', 0.01);
% ds.turbulence_intensity.data  % Access TI data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    arguments
        ds
        options.noise = 0
        options.thresh = 0
        options.field_name = 'turbulence_intensity'
    end
    
    % Validate input dataset structure
    if ~isstruct(ds)
        error('mhkit:dolfyn: Input ds must be a structure');
    end
    
    % Check for required velocity standard deviation field (from ensemble averaging)
    if isfield(ds, 'U_std')
        % Use pre-calculated velocity standard deviation
        vel_std_data = ds.U_std.data;
        dims = ds.U_std.dims;
        coords = ds.U_std.coords;
        
        % Get corresponding velocity magnitude
        if ~isfield(ds, 'U_mag')
            error('mhkit:dolfyn: U_mag field required when using U_std');
        end
        vel_mag_data = ds.U_mag.data;
        
    elseif isfield(ds, 'vel_std')
        % Use velocity component standard deviations
        vel_std_data = ds.vel_std.data;
        dims = ds.vel_std.dims;
        coords = ds.vel_std.coords;
        
        % Store original data for noise correction (before combining components)
        vel_std_original = vel_std_data;
        original_dims = dims;
        
        % Calculate horizontal velocity standard deviation
        dir_idx = find(strcmp(dims, 'dir'));
        if ~isempty(dir_idx)
            if dir_idx == 1
                u_std = squeeze(vel_std_data(1, :));
                v_std = squeeze(vel_std_data(2, :));
            elseif dir_idx == 2
                u_std = squeeze(vel_std_data(:, 1));
                v_std = squeeze(vel_std_data(:, 2));
            else
                u_std = squeeze(vel_std_data(:, :, 1));
                v_std = squeeze(vel_std_data(:, :, 2));
            end
            vel_std_data = sqrt(u_std.^2 + v_std.^2);
            dims = dims(~strcmp(dims, 'dir'));
        end
        
        % Get velocity magnitude
        if isfield(ds, 'U_mag')
            vel_mag_data = ds.U_mag.data;
        else
            error('mhkit:dolfyn: U_mag field required');
        end
        
    else
        error('mhkit:dolfyn: Dataset must contain U_std or vel_std field from ensemble averaging');
    end
    
    % Handle noise parameter and ensure dimension compatibility
    if isstruct(options.noise) && isfield(options.noise, 'data')
        noise_values = options.noise.data;
    else
        noise_values = options.noise;
    end

    % Apply noise correction with proper dimension handling
    if any(noise_values(:) ~= 0)
        % Get data dimensions
        vel_std_size = size(vel_std_data);
        noise_size = size(noise_values);

        % For noise correction, work with original multi-dimensional data if available
        if exist('vel_std_original', 'var') && exist('original_dims', 'var')
            % Work with the original 4D data for proper noise correction
            orig_size = size(vel_std_original);
            
            % Handle dimension expansion for noise values to match original data
            if length(orig_size) == 4 && length(noise_size) == 2
                % vel_std_original: [time, range_singleton, range_bins, components]
                % noise: [range_bins, components]  
                % Need to expand noise to match original dimensions
                noise_expanded = reshape(noise_values, [1, 1, noise_size(1), noise_size(2)]);
                noise_expanded = repmat(noise_expanded, [orig_size(1), orig_size(2), 1, 1]);
                
                % Apply noise correction to original data
                vel_std_original_corrected = sqrt(max(0, vel_std_original.^2 - noise_expanded.^2));
                
                % Now recalculate horizontal standard deviation from corrected data
                dir_idx = find(strcmp(original_dims, 'dir'));
                if ~isempty(dir_idx)
                    if dir_idx == 1
                        u_std_corr = squeeze(vel_std_original_corrected(1, :));
                        v_std_corr = squeeze(vel_std_original_corrected(2, :));
                    elseif dir_idx == 2  
                        u_std_corr = squeeze(vel_std_original_corrected(:, 1));
                        v_std_corr = squeeze(vel_std_original_corrected(:, 2));
                    else
                        u_std_corr = squeeze(vel_std_original_corrected(:, :, 1));
                        v_std_corr = squeeze(vel_std_original_corrected(:, :, 2));
                    end
                    vel_std_corrected = sqrt(u_std_corr.^2 + v_std_corr.^2);
                else
                    vel_std_corrected = vel_std_original_corrected;
                end

            else
                % Fallback to simple correction
                if isscalar(noise_values)
                    noise_expanded = noise_values;
                else
                    noise_expanded = mean(noise_values(:));  % Use mean noise level
                end
                vel_std_corrected = sqrt(max(0, vel_std_data.^2 - noise_expanded.^2));
            end
        else
            % No original multi-dimensional data available
            if isscalar(noise_values)
                noise_expanded = noise_values;
            else
                noise_expanded = mean(noise_values(:));  % Use mean noise level
            end
            vel_std_corrected = sqrt(max(0, vel_std_data.^2 - noise_expanded.^2));
        end
        
        fprintf('  Noise correction applied successfully\\n');
    else
        fprintf('  No noise correction applied (noise = 0)\\n');
        vel_std_corrected = vel_std_data;
    end
    
    % Calculate turbulence intensity
    TI = vel_std_corrected ./ vel_mag_data;
    
    % Apply velocity threshold masking
    TI(vel_mag_data < options.thresh) = NaN;
    
    % Create output structure
    ds_out = ds;
    ds_out.(options.field_name) = struct();
    ds_out.(options.field_name).data = single(TI);
    ds_out.(options.field_name).dims = dims;
    ds_out.(options.field_name).coords = coords;
    ds_out.(options.field_name).units = "1";
    ds_out.(options.field_name).long_name = "Turbulence Intensity";
    ds_out.(options.field_name).description = sprintf(...
        'TI corrected from noise level of %.3f m/s, threshold %.3f m/s', ...
        mean(noise_values(:)), options.thresh);

end

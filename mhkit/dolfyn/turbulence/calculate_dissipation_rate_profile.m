function ds_out = calculate_dissipation_rate_profile(ds_avg, ds_raw, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate turbulent kinetic energy dissipation rate profile from ADCP data
%
% This function processes an entire ADCP velocity profile to calculate
% dissipation rate at each range bin. It automates the workflow of
% extracting velocity data by range, calculating PSD, computing noise levels,
% and applying the Lumley-Terray 1983 method for each bin. Results are
% automatically concatenated into profile-format output fields.
%
% Parameters
% ------------
%   ds_avg: structure
%       Averaged ADCP dataset structure containing velocity magnitude
%   ds_raw: structure
%       Raw ADCP dataset structure containing 5th beam velocity data
%   vel_b5_field: string
%       Name of the 5th beam velocity field in ds_raw
%       Default: 'vel_b5'
%   U_mag_field: string
%       Name of the velocity magnitude field in ds_avg
%       Default: 'U_mag'
%   freq_range: array [1x2]
%       Frequency range for inertial subrange [f_min, f_max] in Hz
%       Default: [0.2, 0.5]
%   noise_field: string
%       Name of existing noise field in ds_avg (optional)
%       If not provided, noise will be calculated for each range bin
%       Default: '' (calculate noise)
%   pct_fN: numeric
%       Percentage of Nyquist frequency for noise calculation
%       Default: 0.9
%   field_name: string
%       Base name for output fields (will create multiple profile fields)
%       Default: 'dissipation_rate'
%   n_fft: numeric
%       FFT length for PSD calculation
%       Default: floor(ds_avg.attrs.n_bin / 2)
%   apply_qc: logical
%       Apply quality control based on turbulence cascade slope
%       Default: true
%   qc_tolerance: numeric
%       Quality control tolerance for cascade slope (percent difference from -5/3)
%       Default: 0.20 (20%)
%
% Returns
% ---------
%   ds_out: structure
%       Input dataset with added profile fields:
%           ds_out.auto_spectra : PSD profile [range x time x freq]
%           ds_out.noise_profile : Noise level profile [range x time]
%           ds_out.dissipation_rate : Dissipation rate profile [range x time]
%           ds_out.qc_slope : Quality control slope values [range x time]
%           ds_out.qc_mask : Quality control mask [range x time]
%
% Notes
% -----
% This function replaces the manual loop approach for processing ADCP profiles.
% It automatically handles range bin extraction, PSD calculation, noise estimation,
% and dissipation rate calculation for the full water column.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    arguments
        ds_avg
        ds_raw
        options.vel_b5_field = 'vel_b5'
        options.U_mag_field = 'U_mag'
        options.freq_range = [0.2, 0.5]
        options.noise_field = ''
        options.pct_fN = 0.9
        options.field_name = 'dissipation_rate'
        options.n_fft = []
        options.apply_qc = true
        options.qc_tolerance = 0.20
    end

    % Validate input structures
    if ~isstruct(ds_avg)
        error('mhkit:dolfyn:calculate_dissipation_rate_profile: ds_avg must be a structure');
    end

    if ~isstruct(ds_raw)
        error('mhkit:dolfyn:calculate_dissipation_rate_profile: ds_raw must be a structure');
    end

    % Check for required fields
    if ~isfield(ds_raw, options.vel_b5_field)
        error('mhkit:dolfyn:calculate_dissipation_rate_profile: ds_raw must contain field: %s', options.vel_b5_field);
    end

    if ~isfield(ds_avg, options.U_mag_field)
        error('mhkit:dolfyn:calculate_dissipation_rate_profile: ds_avg must contain field: %s', options.U_mag_field);
    end

    % Set default n_fft if not provided
    if isempty(options.n_fft)
        if isfield(ds_avg, 'attrs') && isfield(ds_avg.attrs, 'n_bin')
            options.n_fft = floor(ds_avg.attrs.n_bin / 2);
        else
            options.n_fft = 256;  % Default fallback
        end
    end

    % Get range dimensions
    vel_b5_data = ds_raw.(options.vel_b5_field);
    if ~isfield(vel_b5_data, 'coords') || ~isfield(vel_b5_data.coords, 'range_b5')
        error('mhkit:dolfyn:calculate_dissipation_rate_profile: %s must contain range_b5 coordinates', options.vel_b5_field);
    end

    range_coords = vel_b5_data.coords.range_b5;
    n_range = length(range_coords);

    % Initialize storage for results
    spec_profile = cell(n_range, 1);
    noise_profile = cell(n_range, 1);
    dissipation_profile = cell(n_range, 1);

    % Process each range bin
    processed_count = 0;

    for r = 1:n_range
        try
            % Extract velocity data for this range bin
            [vel_b5_r, ~] = dolfyn_select(ds_raw.(options.vel_b5_field), 'range_b5', range_coords(r), 'method', 'nearest');

            % Calculate PSD for this range bin
            ds_temp = power_spectral_density(ds_avg, vel_b5_r, ...
                'freq_units', 'Hz', ...
                'n_fft', options.n_fft, ...
                'field_name', sprintf('auto_spectra_r%d', r));

            % Calculate or use existing noise level
            if isempty(options.noise_field)
                % Calculate noise level for this range bin
                ds_temp = calculate_doppler_noise_level(ds_temp, ...
                    'psd_field', sprintf('auto_spectra_r%d', r), ...
                    'pct_fN', options.pct_fN, ...
                    'field_name', sprintf('noise_r%d', r));
                noise_data = ds_temp.(sprintf('noise_r%d', r)).data;
            else
                % Use existing noise field (extract for this range if needed)
                if isfield(ds_avg, options.noise_field)
                    noise_field_data = ds_avg.(options.noise_field);
                    if isfield(noise_field_data, 'data')
                        noise_data = noise_field_data.data;
                        % If noise is a profile, extract this range
                        if size(noise_data, 1) == n_range
                            noise_data = noise_data(r, :);
                        elseif numel(noise_data) == 1
                            % Scalar noise - use as is
                        else
                            % Use mean if dimensions don't match
                            noise_data = mean(noise_data, 'omitnan');
                        end
                    else
                        noise_data = noise_field_data;
                    end
                else
                    error('mhkit:dolfyn:calculate_dissipation_rate_profile: Specified noise field %s not found', options.noise_field);
                end
            end

            % Extract velocity magnitude for this range bin
            if isfield(ds_avg, options.U_mag_field)
                U_mag_field_data = ds_avg.(options.U_mag_field);

                % Check if U_mag has range dimension
                if isfield(U_mag_field_data, 'coords') && isfield(U_mag_field_data.coords, 'range')
                    [U_mag_r, ~] = dolfyn_select(ds_avg.(options.U_mag_field), 'range', range_coords(r), 'method', 'nearest');
                else
                    % U_mag might be 1D time series - use as is
                    U_mag_r = U_mag_field_data.data;
                end

                % Create temporary 1D U_mag field for this range bin
                u_field_name = sprintf('U_mag_r%d', r);
                ds_temp.(u_field_name) = struct();
                ds_temp.(u_field_name).data = U_mag_r;
                ds_temp.(u_field_name).dims = {'time'};
                ds_temp.(u_field_name).coords = struct();

                % Calculate dissipation rate
                ds_temp = calculate_dissipation_rate_LT83(ds_temp, ...
                    'psd_field', sprintf('auto_spectra_r%d', r), ...
                    'U_mag_field', u_field_name, ...
                    'freq_range', options.freq_range, ...
                    'noise', noise_data, ...
                    'field_name', sprintf('dissipation_r%d', r));

                % Store results
                if isfield(ds_temp, sprintf('auto_spectra_r%d', r))
                    spec_profile{r} = ds_temp.(sprintf('auto_spectra_r%d', r));
                end

                if isempty(options.noise_field) && isfield(ds_temp, sprintf('noise_r%d', r))
                    noise_profile{r} = ds_temp.(sprintf('noise_r%d', r));
                else
                    % Create noise structure for consistency
                    noise_profile{r} = struct();
                    noise_profile{r}.data = noise_data;
                    noise_profile{r}.dims = {'time'};
                    noise_profile{r}.coords = struct();
                end

                if isfield(ds_temp, sprintf('dissipation_r%d', r))
                    dissipation_profile{r} = ds_temp.(sprintf('dissipation_r%d', r));
                end

                processed_count = processed_count + 1;
            else
                fprintf('  Warning: No velocity magnitude data for range bin %d\n', r);
            end

        catch ME
            fprintf('  Warning: Failed to process range bin %d: %s\n', r, ME.message);
            continue;
        end
    end

    ds_out = ds_avg;

    % Concatenate PSD profile
    if ~isempty(spec_profile) && ~isempty(spec_profile{1})
        ds_out = concatenate_profile_results(ds_out, spec_profile, 'auto_spectra', range_coords);
    end

    % Concatenate noise profile
    if ~isempty(noise_profile) && ~isempty(noise_profile{1})
        ds_out = concatenate_profile_results(ds_out, noise_profile, 'noise_profile', range_coords);
    end

    % Concatenate dissipation rate profile
    if ~isempty(dissipation_profile) && ~isempty(dissipation_profile{1})
        ds_out = concatenate_profile_results(ds_out, dissipation_profile, options.field_name, range_coords);
    end

    % Apply quality control if requested
    if options.apply_qc && isfield(ds_out, 'auto_spectra')
        % Calculate cascade slope for quality control
        [slope, ~] = check_turbulence_cascade_slope(ds_out.auto_spectra, 'freq_range', options.freq_range);

        % Check that percent difference from -5/3 is within tolerance
        target_slope = -5.0 / 3.0;
        qc_mask = abs((slope - target_slope) / target_slope) <= options.qc_tolerance;

        % Store QC results
        ds_out.qc_slope = struct();
        ds_out.qc_slope.data = slope;
        ds_out.qc_slope.dims = {'range', 'time'};
        ds_out.qc_slope.coords = struct('range', range_coords);
        ds_out.qc_slope.units = "dimensionless";
        ds_out.qc_slope.long_name = "Turbulence Cascade Slope";
        ds_out.qc_slope.description = sprintf('Spectral slope in frequency range [%.3f, %.3f] Hz', ...
            options.freq_range(1), options.freq_range(2));

        ds_out.qc_mask = struct();
        ds_out.qc_mask.data = qc_mask;
        ds_out.qc_mask.dims = {'range', 'time'};
        ds_out.qc_mask.coords = struct('range', range_coords);
        ds_out.qc_mask.units = "dimensionless";
        ds_out.qc_mask.long_name = "Quality Control Mask";
        ds_out.qc_mask.description = sprintf('QC mask: slope within %.0f%% of -5/3', ...
            options.qc_tolerance * 100);

        % Apply QC mask to dissipation rate
        if isfield(ds_out, options.field_name)
            ds_out.(options.field_name).data(~qc_mask) = NaN;

            % Update metadata to indicate QC was applied
            ds_out.(options.field_name).qc_applied = true;
            ds_out.(options.field_name).qc_tolerance = options.qc_tolerance;
            ds_out.(options.field_name).qc_target_slope = target_slope;
        end

        % QC statistics
        valid_estimates = sum(qc_mask(:));
        total_estimates = numel(qc_mask);

    end

end


function ds_out = concatenate_profile_results(ds_in, profile_cell, field_name, range_coords)
    % Helper function to concatenate cell array results into profile structure

    ds_out = ds_in;

    % Find first non-empty entry to get dimensions
    first_valid = find(~cellfun(@isempty, profile_cell), 1);
    if isempty(first_valid)
        fprintf('Warning: No valid data found for field %s\n', field_name);
        return;
    end

    sample_data = profile_cell{first_valid};

    % Get dimensions from sample
    if isfield(sample_data, 'data')
        sample_size = size(sample_data.data);
        data_dims = length(sample_size);

        % Initialize output array
        if data_dims == 1
            % 1D time series -> [range x time]
            n_range = length(profile_cell);
            n_time = sample_size(1);
            output_data = NaN(n_range, n_time, 'single');
            output_dims = {'range', 'time'};

            % Fill data
            for r = 1:n_range
                if ~isempty(profile_cell{r}) && isfield(profile_cell{r}, 'data')
                    output_data(r, :) = profile_cell{r}.data(:)';
                end
            end

        elseif data_dims == 2
            % 2D [time x freq] -> [range x time x freq]
            n_range = length(profile_cell);
            n_time = sample_size(1);
            n_freq = sample_size(2);
            output_data = NaN(n_range, n_time, n_freq, 'single');
            output_dims = {'range', 'time', 'freq'};

            % Fill data
            for r = 1:n_range
                if ~isempty(profile_cell{r}) && isfield(profile_cell{r}, 'data')
                    output_data(r, :, :) = profile_cell{r}.data;
                end
            end

        else
            error('mhkit:dolfyn:calculate_dissipation_rate_profile: Cannot concatenate %dD data', data_dims);
        end

        % Create output structure
        ds_out.(field_name) = struct();
        ds_out.(field_name).data = output_data;
        ds_out.(field_name).dims = output_dims;

        % Set up coordinates
        ds_out.(field_name).coords = struct();
        ds_out.(field_name).coords.range = range_coords;

        % Copy other coordinates from sample if they exist
        if isfield(sample_data, 'coords')
            coord_fields = fieldnames(sample_data.coords);
            for i = 1:length(coord_fields)
                coord_name = coord_fields{i};
                if ~strcmp(coord_name, 'range')  % Don't overwrite range
                    ds_out.(field_name).coords.(coord_name) = sample_data.coords.(coord_name);
                end
            end
        end

        % Copy metadata from sample
        metadata_fields = {'units', 'long_name', 'standard_name', 'description', 'method'};
        for i = 1:length(metadata_fields)
            if isfield(sample_data, metadata_fields{i})
                ds_out.(field_name).(metadata_fields{i}) = sample_data.(metadata_fields{i});
            end
        end

    else
        fprintf('Warning: Sample data for %s does not contain data field\n', field_name);
    end
end

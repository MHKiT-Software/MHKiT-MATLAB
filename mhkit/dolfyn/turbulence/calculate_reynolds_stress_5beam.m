function ds_out = calculate_reynolds_stress_5beam(ds, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate Reynolds stress using 5-beam variance method for ADCP data.
%
% This function implements the exact algorithms from MHKiT-Python dolfyn module
% for 5-beam Reynolds stress and TKE calculations, following the methods
% from Dewey & Stringer (2007) and Guerra & Thomson (2017).
%
% Parameters
% ------------
%   ds: structure
%       ADCP dataset with 5th beam velocity data
%   vel_field: string, default 'vel'
%       Name of the 4-beam velocity field
%   vel_b5_field: string, default 'vel_b5'  
%       Name of the 5th beam velocity field
%   tke_only: logical, default false
%       If true, only calculate TKE components
%   noise: numeric, default 0
%       Doppler noise level [m/s]
%   beam_angle: numeric, default 20
%       ADCP beam angle in degrees
%   orientation: string, default 'down'
%       Instrument orientation ('up' or 'down')
%   inst_make: string, default 'rdi'
%       Instrument manufacturer ('rdi' or 'nortek')
%   align_with_shear: logical, default true
%       If true, align stress component dimensions with velocity shear
%       (trims range bins to match centered difference output)
%
% Returns
% ---------
%   ds_out: structure
%       Dataset with tke_vec and/or stress_vec fields matching MHKiT-Python output
%
% References
% ----------
% Dewey, R., and S. Stringer. "Reynolds stresses and turbulent kinetic
% energy estimates from various ADCP beam configurations: Theory." J. of
% Phys. Ocean (2007): 1-35.
%
% Guerra, Maricarmen, and Jim Thomson. "Turbulence measurements from
% five-beam acoustic Doppler current profilers." Journal of Atmospheric
% and Oceanic Technology 34.6 (2017): 1267-1284.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    arguments
        ds
        options.vel_field = 'vel'
        options.vel_b5_field = 'vel_b5'
        options.tke_only = false
        options.noise = 0
        options.beam_angle = 20
        options.orientation = 'down'
        options.inst_make = 'rdi'
        options.align_with_shear = true  % Align stress components with velocity shear dimensions
    end
    
    % Basic validation
    if ~isstruct(ds)
        error('mhkit:dolfyn:calculate_reynolds_stress_5beam: Input ds must be a structure');
    end
    
    % Check for required velocity fields
    if ~isfield(ds, options.vel_field)
        error('mhkit:dolfyn:calculate_reynolds_stress_5beam: Dataset must contain 4-beam velocity field: %s', options.vel_field);
    end
    
    if ~isfield(ds, options.vel_b5_field)
        error('mhkit:dolfyn:calculate_reynolds_stress_5beam: Dataset must contain 5th beam velocity field: %s', options.vel_b5_field);
    end
    
    % Validate inputs
    if options.beam_angle <= 0 || options.beam_angle >= 90
        error('mhkit:dolfyn:calculate_reynolds_stress_5beam: Beam angle must be between 0 and 90 degrees (got %.1f°)', options.beam_angle);
    end
    
    if ~ismember(lower(options.orientation), {'up', 'down'})
        error('mhkit:dolfyn:calculate_reynolds_stress_5beam: Orientation must be ''up'' or ''down''');
    end
    
    if ~ismember(lower(options.inst_make), {'rdi', 'nortek'})
        error('mhkit:dolfyn:calculate_reynolds_stress_5beam: inst_make must be ''rdi'' or ''nortek''');
    end
    
    % Get velocity data
    vel_data_raw = ds.(options.vel_field).data;
    vel_b5_data_raw = ds.(options.vel_b5_field).data;
    
    % Determine data dimensions
    vel_size = size(vel_data_raw);
    vel_b5_size = size(vel_b5_data_raw);
    
    % Handle MHKiT-MATLAB data structure format
    % Expected: vel = [time, singleton, range, beam], vel_b5 = [time, singleton, range]
    % Transform to: vel = [range, time, beam], vel_b5 = [range, time]
    
    if length(vel_size) == 4 && vel_size(2) == 1
        % Format: [time, singleton, range, beam] -> [range, time, beam]
        vel_data = squeeze(permute(vel_data_raw, [3, 1, 4, 2]));  % [range, time, beam]
        n_time = vel_size(1);
        n_range = vel_size(3);
        n_beam = vel_size(4);
        
        if n_beam < 4
            error('mhkit:dolfyn:calculate_reynolds_stress_5beam: Expected at least 4 beams, got %d', n_beam);
        end
        
    elseif length(vel_size) == 3
        % Format: [range, time, beam] (already correct)
        vel_data = vel_data_raw;
        [n_range, n_time, n_beam] = size(vel_data);
        
        if n_beam < 4
            error('mhkit:dolfyn:calculate_reynolds_stress_5beam: Expected at least 4 beams, got %d', n_beam);
        end
    else
        error('mhkit:dolfyn:calculate_reynolds_stress_5beam: Unexpected velocity data format. Expected [time, singleton, range, beam] or [range, time, beam], got [%s]', num2str(vel_size));
    end
    
    if length(vel_b5_size) == 3 && vel_b5_size(2) == 1
        % Format: [time, singleton, range] -> [range, time]
        vel_b5_data = squeeze(permute(vel_b5_data_raw, [3, 1, 2]));  % [range, time]
        
        % Verify dimensions match
        if size(vel_b5_data, 1) ~= n_range || size(vel_b5_data, 2) ~= n_time
            error('mhkit:dolfyn:calculate_reynolds_stress_5beam: 5th beam dimensions [%d, %d] do not match velocity data [%d, %d]', ...
                size(vel_b5_data, 1), size(vel_b5_data, 2), n_range, n_time);
        end
        
    elseif length(vel_b5_size) == 2
        % Format: [range, time] (already correct)
        vel_b5_data = vel_b5_data_raw;
        
        if size(vel_b5_data, 1) ~= n_range || size(vel_b5_data, 2) ~= n_time
            error('mhkit:dolfyn:calculate_reynolds_stress_5beam: 5th beam dimensions [%d, %d] do not match velocity data [%d, %d]', ...
                size(vel_b5_data, 1), size(vel_b5_data, 2), n_range, n_time);
        end
    else
        error('mhkit:dolfyn:calculate_reynolds_stress_5beam: Unexpected 5th beam data format. Expected [time, singleton, range] or [range, time], got [%s]', num2str(vel_b5_size));
    end
    
    % Get tilt angles using MHKiT-Python's manufacturer-specific approach
    if isfield(ds, 'pitch') && isfield(ds.pitch, 'data') && ...
       isfield(ds, 'roll') && isfield(ds.roll, 'data')
        
        if strcmpi(options.inst_make, 'rdi')
            % RDI: phi2 = pitch, phi3 = roll
            phi2 = deg2rad(mhkit_nanmean(ds.pitch.data(:)));
            phi3 = deg2rad(mhkit_nanmean(ds.roll.data(:)));
        else
            % Nortek: phi2 = roll, phi3 = -pitch (note the negation!)
            phi2 = deg2rad(mhkit_nanmean(ds.roll.data(:)));
            phi3 = -deg2rad(mhkit_nanmean(ds.pitch.data(:)));  % Critical negation for Nortek
        end
    else
        phi2 = 0;  % No tilt
        phi3 = 0;
    end
    
    % Determine beam order following MHKiT-Python's manufacturer-specific logic
    % MATLAB 1-based: [0, 1, 2, 3] in MHKiT-Python
    if strcmpi(options.inst_make, 'rdi')
        if strcmpi(options.orientation, 'down')
            beam_order = [1, 2, 3, 4];
        else  % up
            beam_order = [1, 2, 4, 3];
        end
    else  % nortek
        if strcmpi(options.orientation, 'down')
            beam_order = [3, 1, 4, 2];
        else  % up
            beam_order = [1, 3, 4, 2];
        end
    end
    
    % Combine 4-beam and 5th beam data: [5, range, time]
    beam_vel = zeros(5, n_range, n_time);
    for i = 1:4
        beam_vel(i, :, :) = squeeze(vel_data(:, :, beam_order(i)));
    end
    beam_vel(5, :, :) = vel_b5_data;  % 5th beam
    
    % Calculate variance along time dimension for each beam and range bin
    bp2_ = zeros(5, n_range);
    for beam = 1:5
        for r = 1:n_range
            bp2_(beam, r) = var(squeeze(beam_vel(beam, r, :)), 'omitnan');
        end
    end
    
    % Apply noise correction
    if options.noise > 0
        bp2_ = bp2_ - options.noise^2;
        bp2_(bp2_ < 0) = NaN;  % Remove negative variances
    end
    
    % 5-beam TKE calculations using exact MHKiT-Python formulas
    th = deg2rad(options.beam_angle);
    denm = -4 * sin(th)^6 * cos(th)^2;
    
    % u'u' component (MHKiT-Python dolfyn/adp/turbulence.py line 688-694)
    upup_ = (-2 * sin(th)^4 * cos(th)^2 .* ...
             (bp2_(2, :) + bp2_(1, :) - 2 * cos(th)^2 .* bp2_(5, :)) ...
             + 2 * sin(th)^5 * cos(th) * phi3 .* (bp2_(2, :) - bp2_(1, :))) ./ denm;
    
    % v'v' component (MHKiT-Python dolfyn/adp/turbulence.py line 696-704)
    vpvp_ = (-2 * sin(th)^4 * cos(th)^2 .* ...
             (bp2_(4, :) + bp2_(3, :) - 2 * cos(th)^2 .* bp2_(5, :)) ...
             - 2 * sin(th)^4 * cos(th)^2 * phi3 .* (bp2_(2, :) - bp2_(1, :)) ...
             + 2 * sin(th)^3 * cos(th)^3 * phi3 .* (bp2_(2, :) - bp2_(1, :)) ...
             - 2 * sin(th)^5 * cos(th) * phi2 .* (bp2_(4, :) - bp2_(3, :))) ./ denm;
    
    % w'w' component (MHKiT-Python dolfyn/adp/turbulence.py line 706-716) - CORRECTED
    wpwp_ = (-2 * sin(th)^5 * cos(th) .* ...
             (bp2_(2, :) - bp2_(1, :) ...
              + 2 * sin(th)^5 * cos(th) * phi2 .* (bp2_(4, :) - bp2_(3, :)) ...
              - 4 * sin(th)^6 * cos(th)^2 .* bp2_(5, :))) ./ denm;
    
    % Dimension Alignment: Match TKE components with velocity shear coordinate system
    if options.align_with_shear
        % When velocity shear is calculated using centered difference, it reduces
        % the range dimension by 2 (loses first and last bins). To ensure compatibility
        % with velocity shear calculations, we align the TKE components to the
        % same range coordinate system that would result from centered difference.
        
        % For centered difference, we lose the first and last range bins
        % So we trim the TKE components to match: [28] -> [26] range bins
        if n_range >= 3  % Ensure we have enough bins for centered difference
            range_start = 2;  % Skip first bin (MATLAB 1-based indexing)
            range_end = n_range - 1;  % Skip last bin
            n_range_shear = range_end - range_start + 1;  % 26 range bins
            
            % Trim TKE components to match shear dimensions
            upup_aligned = upup_(range_start:range_end);
            vpvp_aligned = vpvp_(range_start:range_end);
            wpwp_aligned = wpwp_(range_start:range_end);
            
            % Create aligned range coordinate (matching what dudz/dvdz would produce)
            range_original = ds.(options.vel_field).coords.range;
            range_tke_aligned = range_original(range_start:range_end);
            
            % Stack TKE components with aligned dimensions
            tke_components = zeros(n_range_shear, 1, 3);
            tke_components(:, 1, 1) = upup_aligned;   % upup_ component
            tke_components(:, 1, 2) = vpvp_aligned;   % vpvp_ component
            tke_components(:, 1, 3) = wpwp_aligned;   % wpwp_ component
        else
            error('mhkit:dolfyn:calculate_reynolds_stress_5beam: Insufficient range bins (%d) for TKE dimension alignment', n_range);
        end
    else
        % No alignment - use original TKE components
        tke_components = zeros(n_range, 1, 3);
        tke_components(:, 1, 1) = upup_;   % upup_ component
        tke_components(:, 1, 2) = vpvp_;   % vpvp_ component
        tke_components(:, 1, 3) = wpwp_;   % wpwp_ component
        range_tke_aligned = ds.(options.vel_field).coords.range;
    end
    
    % Create output structure matching expected format
    ds_out = ds;
    ds_out.tke_vec = struct();
    ds_out.tke_vec.data = single(tke_components);
    ds_out.tke_vec.dims = {'range', 'time', 'tke'};
    ds_out.tke_vec.coords = struct();
    ds_out.tke_vec.coords.tke = {'upup_', 'vpvp_', 'wpwp_'};
    ds_out.tke_vec.coords.range = range_tke_aligned;  % Use aligned range coordinate
    ds_out.tke_vec.coords.time = 1;  % Single time point since this is ensemble calculation
    ds_out.tke_vec.units = "m2 s-2";
    ds_out.tke_vec.long_name = "TKE Vector";
    ds_out.tke_vec.standard_name = "specific_turbulent_kinetic_energy_of_sea_water";
    
    % Add Reynolds stress components if not TKE-only mode
    if ~options.tke_only
        % Step 1: Create a temporary dataset structure for rotation
        ds_temp = ds;

        % Step 2: Rotate to instrument coordinates
        % Note: The rotate2 function operates on the dataset's velocity field
        % We need to ensure the dataset is in the correct coordinate system first
        current_coord_sys = ds.coord_sys;

        % try
        % Rotate to instrument coordinates if not already there
        if ~strcmpi(current_coord_sys, 'inst')
            ds_temp = rotate2(ds_temp, 'inst');
        end

        % Step 3: Extract rotated velocity data
        % After rotation, velocity data should be in instrument coordinates
        vel_inst_raw = ds_temp.(options.vel_field).data;

        % Handle data format - convert to [range, time, component] format
        if length(size(vel_inst_raw)) == 4 && size(vel_inst_raw, 2) == 1
            % Format: [time, singleton, range, component] -> [range, time, component]
            vel_inst = squeeze(permute(vel_inst_raw, [3, 1, 4, 2]));
        elseif length(size(vel_inst_raw)) == 3
            % Format: [range, time, component] (already correct)
            vel_inst = vel_inst_raw;
        else
            error('mhkit:dolfyn:calculate_reynolds_stress_5beam: Unexpected rotated velocity data format');
        end

        % Step 4: Extract u and v components (first two components after rotation)
        u_vel_inst = squeeze(vel_inst(:, :, 1));  % u component [range, time]
        v_vel_inst = squeeze(vel_inst(:, :, 2));  % v component [range, time]

        % Step 5: Apply detrending to velocity data (following MHKiT-Python approach)
        u_vel_detrended = zeros(size(u_vel_inst));
        v_vel_detrended = zeros(size(v_vel_inst));

        for r = 1:n_range
            % Detrend each range bin's time series using mhkit_detrend_array
            u_vel_detrended(r, :) = mhkit_detrend_array(u_vel_inst(r, :));
            v_vel_detrended(r, :) = mhkit_detrend_array(v_vel_inst(r, :));
        end

        % Step 6: Calculate u'v' covariance from detrended data
        upvp_ = zeros(n_range, 1);
        for r = 1:n_range
            % Calculate covariance: E[u' * v'] = mean(u' * v')
            upvp_(r) = mean(u_vel_detrended(r, :) .* v_vel_detrended(r, :), 'omitnan');
        end
        
        % u'w' and v'w' components using MHKiT-Python formulas (lines 743-755)
        % Note: bp2_ is [5 x range], we compute as row vectors then transpose to match upvp_
        upwp_ = (sin(th)^5 * cos(th) .* (bp2_(2, :) - bp2_(1, :)) ...
                 + 2 * sin(th)^4 * cos(th) * 2 * phi3 .* (bp2_(2, :) + bp2_(1, :)) ...
                 - 4 * sin(th)^4 * cos(th) * 2 * phi3 .* bp2_(5, :) ...
                 - 4 * sin(th)^6 * cos(th) * 2 * phi2 .* upvp_') ./ denm;
        
        vpwp_ = (sin(th)^5 * cos(th) .* (bp2_(4, :) - bp2_(3, :)) ...
                 - 2 * sin(th)^4 * cos(th) * 2 * phi2 .* (bp2_(4, :) + bp2_(3, :)) ...
                 + 4 * sin(th)^4 * cos(th) * 2 * phi2 .* bp2_(5, :) ...
                 + 4 * sin(th)^6 * cos(th) * 2 * phi3 .* upvp_') ./ denm;
        
        % Transpose to match upvp_ dimensions [n_range x 1]
        upwp_ = upwp_';
        vpwp_ = vpwp_';
        
        % DIMENSION ALIGNMENT: Match velocity shear coordinate system
        if options.align_with_shear
            % When velocity shear is calculated using centered difference, it reduces
            % the range dimension by 2 (loses first and last bins). To ensure compatibility
            % with velocity shear calculations, we align the stress components to the
            % same range coordinate system that would result from centered difference.

            % For centered difference, we lose the first and last range bins
            % So we trim the stress components to match: [28] -> [26] range bins
            if n_range >= 3  % Ensure we have enough bins for centered difference
                range_start = 2;  % Skip first bin (MATLAB 1-based indexing)
                range_end = n_range - 1;  % Skip last bin
                n_range_shear = range_end - range_start + 1;  % 26 range bins
                
                % Trim stress components to match shear dimensions
                upvp_aligned = upvp_(range_start:range_end);
                upwp_aligned = upwp_(range_start:range_end);
                vpwp_aligned = vpwp_(range_start:range_end);
                
                % Create aligned range coordinate (matching what dudz/dvdz would produce)
                range_original = ds.(options.vel_field).coords.range;
                range_aligned = range_original(range_start:range_end);
                
                % Stack stress components with aligned dimensions
                stress_components = zeros(n_range_shear, 1, 3);
                stress_components(:, 1, 1) = upvp_aligned;  % upvp_ component
                stress_components(:, 1, 2) = upwp_aligned;  % upwp_ component  
                stress_components(:, 1, 3) = vpwp_aligned;  % vpwp_ component
            else
                error('mhkit:dolfyn:calculate_reynolds_stress_5beam: Insufficient range bins (%d) for dimension alignment', n_range);
            end
        else
            % No alignment - use original stress components
            stress_components = zeros(n_range, 1, 3);
            stress_components(:, 1, 1) = upvp_;  % upvp_ component
            stress_components(:, 1, 2) = upwp_;  % upwp_ component  
            stress_components(:, 1, 3) = vpwp_;  % vpwp_ component
            range_aligned = ds.(options.vel_field).coords.range;
        end
        
        ds_out.stress_vec = struct();
        ds_out.stress_vec.data = single(stress_components);
        ds_out.stress_vec.dims = {'range', 'time', 'tau'};
        ds_out.stress_vec.coords = struct();
        ds_out.stress_vec.coords.tau = {'upvp_', 'upwp_', 'vpwp_'};
        ds_out.stress_vec.coords.range = range_aligned;  % Use aligned range coordinate
        ds_out.stress_vec.coords.time = 1;  % Single time point since this is ensemble calculation
        ds_out.stress_vec.units = "m2 s-2";
        ds_out.stress_vec.long_name = "Specific Reynolds Stress Vector";
        ds_out.stress_vec.standard_name = "specific_momentum_flux_in_sea_water";
        if options.align_with_shear
            ds_out.stress_vec.description = sprintf(...
                '5-beam Reynolds stress using %s %s-facing configuration with %.1f° beam angle and %.3f° tilt (shear-aligned)', ...
                options.inst_make, options.orientation, options.beam_angle, rad2deg(sqrt(phi2^2 + phi3^2)));
            ds_out.stress_vec.shear_alignment = 'centered';  % Flag indicating alignment with centered difference
        else
            ds_out.stress_vec.description = sprintf(...
                '5-beam Reynolds stress using %s %s-facing configuration with %.1f° beam angle and %.3f° tilt', ...
                options.inst_make, options.orientation, options.beam_angle, rad2deg(sqrt(phi2^2 + phi3^2)));
            ds_out.stress_vec.shear_alignment = 'none';  % Flag indicating no alignment
        end
        
    end
    
    % Add comprehensive metadata
    if options.align_with_shear
        ds_out.tke_vec.description = sprintf(...
            '5-beam TKE using %s %s-facing configuration with %.1f° beam angle and %.3f° tilt (shear-aligned)', ...
            options.inst_make, options.orientation, options.beam_angle, rad2deg(sqrt(phi2^2 + phi3^2)));
        ds_out.tke_vec.shear_alignment = 'centered';  % Flag indicating alignment with centered difference
    else
        ds_out.tke_vec.description = sprintf(...
            '5-beam TKE using %s %s-facing configuration with %.1f° beam angle and %.3f° tilt', ...
            options.inst_make, options.orientation, options.beam_angle, rad2deg(sqrt(phi2^2 + phi3^2)));
        ds_out.tke_vec.shear_alignment = 'none';  % Flag indicating no alignment
    end
    ds_out.tke_vec.method = '5-beam variance method (Dewey & Stringer 2007; Guerra & Thomson 2017)';
    ds_out.tke_vec.beam_angle = options.beam_angle;
    ds_out.tke_vec.tilt_angles = [rad2deg(phi2), rad2deg(phi3)];
    ds_out.tke_vec.noise_correction = options.noise;
    ds_out.tke_vec.beam_order = beam_order;
end

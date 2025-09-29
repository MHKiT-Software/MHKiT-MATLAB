classdef Dolfyn_Test_Analysis_Workflow < matlab.unittest.TestCase

    properties
        test_data_path
        tolerance
        adv_data
        adcp_data
        tidal_data
        tidal_processed
    end

    methods (TestMethodSetup)
        function setupTestData(testCase)
            testCase.tolerance = 1e-6;
            testCase.test_data_path = '../../examples/data/dolfyn';

            adv_file = fullfile(testCase.test_data_path, 'vector_data01.VEC');
            verifyTrue(testCase, exist(adv_file, 'file') == 2, ...
                sprintf('ADV test data file must exist: %s', adv_file));
            testCase.adv_data = dolfyn_read(adv_file);

            adcp_file = fullfile(testCase.test_data_path, 'BenchFile01.ad2cp');
            verifyTrue(testCase, exist(adcp_file, 'file') == 2, ...
                sprintf('ADCP test data file must exist: %s', adcp_file));
            testCase.adcp_data = dolfyn_read(adcp_file);

            tidal_file = fullfile(testCase.test_data_path, 'Sig1000_tidal.ad2cp');
            verifyTrue(testCase, exist(tidal_file, 'file') == 2, ...
                sprintf('Tidal ADCP test data file must exist: %s', tidal_file));
            testCase.tidal_data = dolfyn_read(tidal_file, 'nens', [10000, 19000]);

            testCase.tidal_processed = testCase.processForTurbulence(testCase.tidal_data);
        end
    end

    methods
        function ds_processed = processForTurbulence(~, ds)
            % Execute the same steps as adcp_example.m
            ds = set_range_offset(ds, 0.6);
            ds = water_depth_from_pressure(ds, 'salinity', 31);
            ds = remove_surface_interference(ds);
            ds = correlation_filter(ds, 'thresh', 50);
            ds = set_declination(ds, 15.8);
            ds = rotate2(ds, 'earth');

            sampling_frequency_hz = ds.attrs.fs;
            averaging_period_seconds = 300;  % Smaller period to work with 100 samples
            averaging_samples = int32(round(averaging_period_seconds * sampling_frequency_hz));
            ds_processed = average_by_dimension(ds, averaging_samples, 'time');
            ds_processed = calculate_horizontal_speed_and_direction(ds_processed);
        end
    end

    methods (Test)

function test_turbulence_intensity_basic(testCase)
    result = calculate_turbulence_intensity(testCase.tidal_processed, 'field_name', 'TI');
    verifyTrue(testCase, isstruct(result), 'Result should be a structure');
    verifyTrue(testCase, isfield(result, 'TI'), 'Result should contain TI field');

    if isfield(result, 'TI') && isfield(result.TI, 'data')
        ti_data = result.TI.data;
        verifyTrue(testCase, isvector(ti_data) || ismatrix(ti_data), 'TI data should be vector or matrix');
        verifyGreaterThan(testCase, numel(ti_data), 0, 'TI should contain data');
        verifyTrue(testCase, isnumeric(ti_data), 'TI data should be numeric');
    end
end

function test_power_spectral_density_workflow(testCase)
    rng = 5;
    [vel_up, ~] = dolfyn_select(testCase.tidal_data.vel_b5, 'range_b5', rng, 'method', 'nearest');
    [U_data, ~] = dolfyn_select(testCase.tidal_processed.U_mag, 'range', rng, 'method', 'nearest');

    result = power_spectral_density(testCase.tidal_processed, vel_up, ...
        'freq_units', 'Hz', ...
        'n_fft', floor(testCase.tidal_processed.attrs.n_bin / 2), ...
        'field_name', 'auto_spectra_5m');

    verifyTrue(testCase, isstruct(result), 'Result should be a structure');
    verifyTrue(testCase, isfield(result, 'auto_spectra_5m'), 'Result should contain PSD field');

    if isfield(result, 'auto_spectra_5m') && isfield(result.auto_spectra_5m, 'data')
        psd_data = result.auto_spectra_5m.data;
        verifyTrue(testCase, all(psd_data(:) >= 0), 'PSD should be non-negative');
        verifyTrue(testCase, all(isfinite(psd_data(:))), 'PSD should be finite');
    end
end

function test_doppler_noise_level_workflow(testCase)
    rng = 5;
    [vel_up, ~] = dolfyn_select(testCase.tidal_data.vel_b5, 'range_b5', rng, 'method', 'nearest');

    ds_with_psd = power_spectral_density(testCase.tidal_processed, vel_up, ...
        'freq_units', 'Hz', ...
        'n_fft', floor(testCase.tidal_processed.attrs.n_bin / 2), ...
        'field_name', 'auto_spectra_5m');

    result = calculate_doppler_noise_level(ds_with_psd, ...
        'psd_field', 'auto_spectra_5m', ...
        'pct_fN', 0.9, ...
        'field_name', 'noise_5m');

    verifyTrue(testCase, isstruct(result), 'Result should be a structure');
    verifyTrue(testCase, isfield(result, 'noise_5m'), 'Result should contain noise field');

    if isfield(result, 'noise_5m') && isfield(result.noise_5m, 'data')
        noise_data = result.noise_5m.data;
        verifyTrue(testCase, all(noise_data(:) >= 0), 'Noise should be non-negative');
        verifyTrue(testCase, all(isfinite(noise_data(:))), 'Noise should be finite');
    end
end

function test_dissipation_rate_workflow(testCase)
    rng = 5;
    [vel_up, ~] = dolfyn_select(testCase.tidal_data.vel_b5, 'range_b5', rng, 'method', 'nearest');
    [U_data, ~] = dolfyn_select(testCase.tidal_processed.U_mag, 'range', rng, 'method', 'nearest');

    ds_with_psd = power_spectral_density(testCase.tidal_processed, vel_up, ...
        'freq_units', 'Hz', ...
        'n_fft', floor(testCase.tidal_processed.attrs.n_bin / 2), ...
        'field_name', 'auto_spectra_5m');

    ds_with_noise = calculate_doppler_noise_level(ds_with_psd, ...
        'psd_field', 'auto_spectra_5m', ...
        'pct_fN', 0.9, ...
        'field_name', 'noise_5m');

    ds_with_noise.U_mag_5m = struct();
    ds_with_noise.U_mag_5m.data = U_data;
    ds_with_noise.U_mag_5m.dims = {'time'};
    ds_with_noise.U_mag_5m.coords = struct();

    f_rng = [0.2, 0.5];
    result = calculate_dissipation_rate_LT83(ds_with_noise, ...
        'psd_field', 'auto_spectra_5m', ...
        'U_mag_field', 'U_mag_5m', ...
        'freq_range', f_rng, ...
        'noise', ds_with_noise.noise_5m, ...
        'field_name', 'dissipation_rate_5m');

    verifyTrue(testCase, isstruct(result), 'Result should be a structure');
    verifyTrue(testCase, isfield(result, 'dissipation_rate_5m'), 'Result should contain dissipation field');
end

function test_reynolds_stress_5beam_workflow(testCase)
    ds_beam = rotate2(testCase.tidal_data, 'beam');

    rng = 5;
    [vel_up, ~] = dolfyn_select(testCase.tidal_data.vel_b5, 'range_b5', rng, 'method', 'nearest');
    ds_with_psd = power_spectral_density(testCase.tidal_processed, vel_up, ...
        'freq_units', 'Hz', ...
        'n_fft', floor(testCase.tidal_processed.attrs.n_bin / 2), ...
        'field_name', 'auto_spectra_5m');
    ds_with_noise = calculate_doppler_noise_level(ds_with_psd, ...
        'psd_field', 'auto_spectra_5m', ...
        'pct_fN', 0.9, ...
        'field_name', 'noise_5m');

    noise_mean = mean(ds_with_noise.noise_5m.data(:), 'omitnan');
    result = calculate_reynolds_stress_5beam(ds_beam, ...
        'noise', noise_mean, ...
        'orientation', 'up', ...
        'beam_angle', 25, ...
        'align_with_shear', true);

    verifyTrue(testCase, isstruct(result), 'Result should be a structure');

    expected_fields = {'tke_vec', 'stress_vec'};
    for i = 1:length(expected_fields)
        verifyTrue(testCase, isfield(result, expected_fields{i}), ...
            sprintf('Result should contain %s field', expected_fields{i}));
    end
end

function test_reynolds_stress_4beam_workflow(testCase)
    % Squeeze processed data to handle dimension mismatches
    ds_squeezed = testCase.tidal_processed;
    if isfield(ds_squeezed, 'vel') && ndims(ds_squeezed.vel.data) > 3
        ds_squeezed.vel.data = squeeze(ds_squeezed.vel.data);
    end

    result = calculate_reynolds_stress_4beam(ds_squeezed);
    verifyTrue(testCase, isstruct(result), 'Result should be a structure');

    result_fields = fieldnames(result);
    stress_found = false;
    for i = 1:length(result_fields)
        if contains(result_fields{i}, 'stress') || contains(result_fields{i}, 'reynolds')
            stress_found = true;
            break;
        end
    end
    verifyTrue(testCase, stress_found, 'Should find stress field in 4-beam result');
end

function test_velocity_shear_workflow(testCase)
    % Squeeze processed data to handle dimension mismatches
    ds_squeezed = testCase.tidal_processed;
    if isfield(ds_squeezed, 'vel') && ndims(ds_squeezed.vel.data) > 3
        ds_squeezed.vel.data = squeeze(ds_squeezed.vel.data);
    end

    result = calculate_velocity_shear(ds_squeezed, ...
        'vel_field', 'vel', ...
        'diff_style', 'centered');

    verifyTrue(testCase, isstruct(result), 'Result should be a structure');

    expected_fields = {'dudz', 'dvdz', 'dwdz'};
    for i = 1:length(expected_fields)
        verifyTrue(testCase, isfield(result, expected_fields{i}), ...
            sprintf('Result should contain %s field', expected_fields{i}));
    end
end

function test_dissipation_rate_profile_workflow(testCase)
    f_rng = [0.2, 0.5];
    result = calculate_dissipation_rate_profile(testCase.tidal_processed, testCase.tidal_data, 'freq_range', f_rng);

    verifyTrue(testCase, isstruct(result), 'Result should be a structure');
    verifyTrue(testCase, isfield(result, 'dissipation_rate'), 'Result should contain dissipation_rate field');
    verifyTrue(testCase, isfield(result, 'qc_mask'), 'Result should contain qc_mask field');
    verifyTrue(testCase, isfield(result, 'qc_slope'), 'Result should contain qc_slope field');
end

function test_data_structure_processed(testCase)
    verifyTrue(testCase, isstruct(testCase.tidal_processed), 'Processed data should be a structure');
    verifyTrue(testCase, isfield(testCase.tidal_processed, 'vel'), 'Should have vel field');
    verifyTrue(testCase, isfield(testCase.tidal_processed, 'U_mag'), 'Should have U_mag field');
    verifyTrue(testCase, isfield(testCase.tidal_processed, 'U_dir'), 'Should have U_dir field');

    if isfield(testCase.tidal_processed, 'vel') && isfield(testCase.tidal_processed.vel, 'data')
        vel_data = testCase.tidal_processed.vel.data;
        vel_size = size(vel_data);
        verifyGreaterThanOrEqual(testCase, length(vel_size), 3, 'Processed velocity should be at least 3D');
        verifyGreaterThan(testCase, vel_size(1), 5, 'Should have time samples');
    end
end

function test_preprocessing_chain(testCase)
    ds = testCase.tidal_data;
    ds = set_range_offset(ds, 0.6);
    verifyTrue(testCase, isstruct(ds), 'set_range_offset should return structure');

    ds = water_depth_from_pressure(ds, 'salinity', 31);
    verifyTrue(testCase, isstruct(ds), 'water_depth_from_pressure should return structure');
    verifyTrue(testCase, isfield(ds, 'depth'), 'Should add depth field');

    ds = remove_surface_interference(ds);
    verifyTrue(testCase, isstruct(ds), 'remove_surface_interference should return structure');

    ds = correlation_filter(ds, 'thresh', 50);
    verifyTrue(testCase, isstruct(ds), 'correlation_filter should return structure');

    ds = set_declination(ds, 15.8);
    verifyTrue(testCase, isstruct(ds), 'set_declination should return structure');

    ds = rotate2(ds, 'earth');
    verifyTrue(testCase, isstruct(ds), 'rotate2 should return structure');

    sampling_frequency_hz = ds.attrs.fs;
    averaging_period_seconds = 300;
    averaging_samples = int32(round(averaging_period_seconds * sampling_frequency_hz));
    ds_avg = average_by_dimension(ds, averaging_samples, 'time');
    verifyTrue(testCase, isstruct(ds_avg), 'average_by_dimension should return structure');

    ds_avg = calculate_horizontal_speed_and_direction(ds_avg);
    verifyTrue(testCase, isstruct(ds_avg), 'calculate_horizontal_speed_and_direction should return structure');
    verifyTrue(testCase, isfield(ds_avg, 'U_mag'), 'Should create U_mag field');
    verifyTrue(testCase, isfield(ds_avg, 'U_dir'), 'Should create U_dir field');
end

function test_input_validation(testCase)
    empty_ds = struct();
    verifyError(testCase, @() calculate_turbulence_intensity(empty_ds), '', ...
        'Should error on empty dataset');

    bad_ds = struct();
    bad_ds.vel = struct();
    verifyError(testCase, @() calculate_turbulence_intensity(bad_ds), '', ...
        'Should error on malformed dataset');
end

    end

end

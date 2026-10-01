function out = convert_to_custom_bands(spsd, bands_per_division, base, use_fft_res_at_bottom)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Convert sound pressure spectral density to custom band spacing based on specified
% parameters.
%
% Parameters
% ------------
%   spsd: struct
%       Sound pressure spectral density structure containing:
%       spsd.data : Spectral density data [Pa^2/Hz or V^2/Hz]
%       spsd.freq : Frequency vector [Hz]
%       spsd.time : Time vector
%   bands_per_division: int
%       The number of bands to divide the spectrum into per increase by a factor
%       of 'base'. A base of 2 and bands_per_division of 3 results in third
%       octaves base 2. Base 10 and bands_per_division of 1000 results in
%       millidecades.
%   base: int
%       Base for the band levels, generally 10 or 2.
%   use_fft_res_at_bottom: logical (optional)
%       In some cases, like millidecades, we do not want to have logarithmically
%       spaced frequency bands across the full spectrum, instead we have the option
%       to have bands that are equal to the FFT bin size. The switch to log spacing
%       is made at the band that has a bandwidth greater than the FFT bin size and
%       such that the frequency space between band center frequencies is at least
%       the FFT bin size. Default is false.
%
% Returns
% ---------
%   out: struct
%       Frequency band-averaged spectral density structure containing:
%       out.data : Band-averaged spectral density data [Units^2/Hz]
%       out.freq : Center frequencies of bands [Hz]
%       out.time : Time vector
%       out.name : Descriptive string
%       out.units : Units string
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    spsd struct
    bands_per_division {mustBeInteger, mustBePositive}
    base {mustBeInteger, mustBePositive}
    use_fft_res_at_bottom logical = false
end

arguments (Output)
    out struct
end

% Validate spsd structure
validate_spsd_struct(spsd, 'convert_to_custom_bands', ...
    'required_fields', {{'data', 'freq', 'time'}});

% Get bands table [lower, center, upper]
bands_table = get_band_table(spsd.freq, bands_per_division, base, use_fft_res_at_bottom);

% Calculate bin indices and fractional weights
[full_pts, partial_pts, weights] = band_power_spectral_density_v3(spsd.freq, bands_table);

% Integrate and average the spectral density
band_spsd_data = band_mean_power_spectral_density_v2(spsd.data, spsd.freq, bands_table, full_pts, partial_pts, weights);

% Build returned structure
out = spsd;
out.data = band_spsd_data;
out.freq = bands_table(:, 2);

end

function bands = get_band_table(freq, bands_per_division, base, use_fft_res_at_bottom)
    fft_bin_size = freq(2) - freq(1);
    bin1_center_freq = freq(1);
    max_freq = freq(end);

    % Generate the log-spaced bands
    [~, band] = create_frequency_bands(bands_per_division, base, freq(1), freq(end));

    band_count_m = 1;
    linear_bin_count = 0;
    max_linear_bin_hz = 0;

    if use_fft_res_at_bottom
        % Find the first band where the bandwidth is greater than or equal to the FFT bin size
        % and the frequency space between band center frequencies is at least the FFT bin size
        diffs = diff(band.center_freq);
        idx_m = find(diffs >= fft_bin_size);
        if length(idx_m) < 2
            band_count_m = 1;
        else
            band_count_m = idx_m(2);
        end
        center_freq = band.center_freq(band_count_m);

        % Now keep counting until the difference between the log spaced
        % center frequency and new frequency is greater than .025
        max_linear_bin_hz = ceil(center_freq / fft_bin_size);
        while (max_linear_bin_hz * fft_bin_size - center_freq) > 0
            band_count_m = band_count_m + 1;
            max_linear_bin_hz = max_linear_bin_hz + 1;
            k_power = band_count_m - 1;
            center_freq = bin1_center_freq * base^(k_power / bands_per_division);
        end

        if (fft_bin_size * max_linear_bin_hz) > max_freq
            max_linear_bin_hz = max_freq / fft_bin_size + 1;
        end

        linear_bin_count = floor(max_linear_bin_hz - 1);
    end

    % Log spaced bands — exclude centers at or beyond max_freq to match Python's
    % exclusive upper-bound behaviour in np.arange inside create_frequency_bands.
    idx = find(band.center_freq > max_linear_bin_hz * fft_bin_size & band.center_freq < max_freq);
    log_bin_count = length(idx);

    % Generate the linear and log frequencies
    total_bins = linear_bin_count + log_bin_count;
    bands = zeros(total_bins, 3);

    if linear_bin_count > 0
        % Linear centers, lower limits, upper limits
        centers = (0:linear_bin_count-1) * fft_bin_size + bin1_center_freq;
        lowers = centers - fft_bin_size / 2;
        uppers = centers + fft_bin_size / 2;

        bands(1:linear_bin_count, 1) = lowers';
        bands(1:linear_bin_count, 2) = centers';
        bands(1:linear_bin_count, 3) = uppers';
    end

    bands(linear_bin_count+1:end, 1) = band.lower_limit(idx);
    bands(linear_bin_count+1:end, 2) = band.center_freq(idx);
    bands(linear_bin_count+1:end, 3) = band.upper_limit(idx);

    if log_bin_count > 0
        bands(end, 3) = max_freq;
    end
end

function [full_pts, partial_pts, weights] = band_power_spectral_density_v3(freq_fft, freq_table)
    fft_bin_size = freq_fft(2) - freq_fft(1);
    num_bands = size(freq_table, 1);

    full_pts = cell(num_bands, 1);
    partial_pts = cell(num_bands, 1);
    weights = cell(num_bands, 1);

    for j = 1:num_bands
        f_lo = freq_table(j, 1);
        f_hi = freq_table(j, 3);

        % FFT bins whose extent falls entirely within the frequency band
        full_mask = ((freq_fft - fft_bin_size / 2) >= f_lo) & ...
                    ((freq_fft + fft_bin_size / 2) <= f_hi);
        full_pts{j} = find(full_mask);

        % FFT bins that overlap the frequency band at all (full or partial)
        overlap_mask = ((freq_fft >= f_lo) & (freq_fft <= f_hi)) | ...
                       (((freq_fft + fft_bin_size / 2) > f_lo) & ...
                        ((freq_fft - fft_bin_size / 2) < f_hi));

        % Partial bins = overlapping minus fully-contained
        partial_pts{j} = setdiff(find(overlap_mask), full_pts{j});
    end

    % Compute fractional weights for partial bins
    for j = 1:num_bands
        num_partial = length(partial_pts{j});
        w = zeros(num_partial, 1);

        for k = 1:num_partial
            idx = partial_pts{j}(k);
            bin_lo = freq_fft(idx) - fft_bin_size / 2;
            bin_hi = freq_fft(idx) + fft_bin_size / 2;

            if bin_lo < freq_table(j, 1)
                % Bin extends below the band's lower edge
                w(k) = (fft_bin_size - (freq_table(j, 1) - bin_lo)) / fft_bin_size;
            elseif bin_hi > freq_table(j, 3)
                % Bin extends above the band's upper edge
                w(k) = (fft_bin_size - (bin_hi - freq_table(j, 3))) / fft_bin_size;
            end
        end
        weights{j} = w;
    end
end

function out_spsd = band_mean_power_spectral_density_v2(input_spsd, freq_fft, freq_table, full_pts, partial_pts, weights)
    fft_bin_size = freq_fft(2) - freq_fft(1);
    nTime = size(input_spsd, 2);
    nBands = size(freq_table, 1);
    out_spsd = zeros(nBands, nTime);

    for j = 1:nBands
        % Contribution from fully-contained FFT bins
        if ~isempty(full_pts{j})
            out_spsd(j, :) = sum(input_spsd(full_pts{j}, :), 1) * fft_bin_size;
        end

        % Contribution from partial FFT bins
        if ~isempty(partial_pts{j})
            weighted_partial = input_spsd(partial_pts{j}, :) .* weights{j};
            out_spsd(j, :) = out_spsd(j, :) + sum(weighted_partial, 1) * fft_bin_size;
        end
    end

    % Take means by dividing by band widths
    band_widths = freq_table(:, 3) - freq_table(:, 1);
    out_spsd = out_spsd ./ band_widths;
end

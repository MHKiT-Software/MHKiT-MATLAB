function df = typecast_spectra_to_mhkit_python(S)
    % Initialize empty python list for column names
    column_names = py.list();

    % Validate frequency is 1D and make row-wise
    freq_dims = size(S.frequency);
    if length(freq_dims) > 2 || (freq_dims(1) ~= 1 && freq_dims(2) ~= 1)
        error('Frequency must be a 1D array. Current dimensions are %s', mat2str(freq_dims));
    end
    freq_length = length(S.frequency);
    if freq_dims(1) ~= 1
        S.frequency = S.frequency(:)';
    end

    % Get initial spectrum dimensions
    [spec_rows, spec_cols] = size(S.spectrum);

    % Case 1: Square matrix matching frequency length
    if spec_rows == spec_cols && spec_rows == freq_length
        warning('Input spectrum is square (%d×%d) and matches frequency length. Unable to determine orientation. Assuming input orientation is correct.', spec_rows, spec_cols);

    % Case 2: Already correct shape (rows match frequency length)
    elseif spec_cols == freq_length
        % No action needed, already correct

    % Case 3: Wrong orientation but can be fixed
    elseif spec_rows == freq_length
        S.spectrum = S.spectrum';
        [spec_rows, spec_cols] = size(S.spectrum);

    % Case 4: Invalid dimensions
    else
        error('Spectrum dimensions (%d×%d) cannot be reconciled with frequency length (%d). Spectrum should be either 1×%d or M×%d.', ...
              spec_rows, spec_cols, freq_length, freq_length, freq_length);
    end

    % Create the frequency index
    freq_index = py.numpy.array(S.frequency);

    % Handle single or multiple spectra
    if spec_rows == 1
        % Single spectrum case
        spectrum_data = py.numpy.array(S.spectrum);
        % Use type as column name, default to 'spectrum' if not provided
        if isfield(S, 'type')
            column_names.append(py.str(S.type));
        else
            column_names.append(py.str('spectrum'));
        end
    else
        % Multiple spectra case
        % Convert all columns to numpy arrays and store column names
        spectrum_data = py.numpy.zeros([spec_rows, spec_cols]);
        for i = 1:spec_rows
            spectrum_data(:,i) = py.numpy.array(S.spectrum(i,:));
            % Generate column name for each spectrum
            if isfield(S, 'type') && length(S.type) >= i
                column_names.append(py.str(S.type{i}));
            else
                column_names.append(py.str(sprintf('spectrum_%d', i)));
            end
        end
    end

    % Create pandas DataFrame with frequency as index
    df = py.pandas.DataFrame(data=spectrum_data, ...
                           index=freq_index, ...
                           columns=column_names);
end

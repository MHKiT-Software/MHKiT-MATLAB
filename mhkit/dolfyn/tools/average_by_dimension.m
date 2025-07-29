function ds_out = average_by_dimension(ds, n_samples, dim_to_find)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Averages data along a specified dimension in a dataset structure
%
% Parameters
% ------------
%     ds: structure
%         Dataset structure containing fields with:
%           .data: N-dimensional array of data
%           .dims: cell array of dimension names
%           .coords (optional): structure of coordinate data
%
%     n_samples: integer
%         Number of samples to average together. Must be positive.
%         For uneven divisions, the function will use floor(dimension_length/n_samples)
%         complete bins and discard remaining samples.
%
%     dim_to_find: string (optional, default='time')
%         Name of the dimension along which to perform averaging.
%         Must exist in at least one field's .dims array.
%
% Returns
% ---------
%     ds_out: structure
%         Output dataset with same structure as input but averaged data:
%           - Data arrays will have reduced size along averaged dimension
%           - Coordinate arrays will be subsampled to match new data size
%           - All other fields and dimensions remain unchanged
%
% Examples
% ---------
%     % Average every 2 time points
%     ds_averaged = average_by_dimension(ds, 2, 'time')
%
%     % Average every 5 points along spatial dimension
%     ds_averaged = average_by_dimension(ds, 5, 'space')
%
% Notes
% -----
%     - If n_samples doesn't divide evenly into the dimension length,
%       the function will issue a warning and discard remaining samples
%     - NaN values are handled appropriately in the averaging
%     - Coordinates are subsampled rather than averaged
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   if nargin < 3
        dim_to_find = 'time';
    end
    if nargin < 2
        error('n_samples is required');
    end

    ds_out = ds;  % Start with a copy of the input

    % Validate dimension name exists in at least one field
    dim_exists = false;
    fields = fieldnames(ds);
    for i = 1:length(fields)
        current_field = fields{i};
        if isfield(ds.(current_field), 'dims') && iscell(ds.(current_field).dims)
            if any(strcmp(ds.(current_field).dims, dim_to_find))
                dim_exists = true;
                break;
            end
        end
    end

    if ~dim_exists
        error('MATLAB:error', 'Dimension "%s" not found in any field', dim_to_find);
    end

    % Process each field that has dimensions
    for i = 1:length(fields)
        current_field = fields{i};
        if isfield(ds.(current_field), 'dims') && iscell(ds.(current_field).dims)
            % Find which dimension matches dim_to_find
            dim_idx = find(strcmp(ds.(current_field).dims, dim_to_find));
            if ~isempty(dim_idx)
                % Get size of the dimension we're averaging
                sz = size(ds.(current_field).data);
                dim_length = sz(dim_idx);

                % Calculate number of complete bins
                n_bins = floor(dim_length / n_samples);

                % Check if n_samples divides evenly
                remainder = mod(dim_length, n_samples);
                if remainder ~= 0
                    warning(['Number of samples (' num2str(n_samples) ') does not divide evenly into length of ' ...
                            current_field ' dimension ' dim_to_find ' (' num2str(dim_length) ...
                            '). Last ' num2str(remainder) ' samples will be discarded.']);
                end

                % Only process if we have at least one complete bin
                if n_bins > 0
                    % Prepare permutation order to move target dimension to first position
                    perm_order = 1:ndims(ds.(current_field).data);
                    perm_order(1) = dim_idx;
                    perm_order(dim_idx) = 1;

                    % Permute data to put target dimension first
                    data = permute(ds.(current_field).data, perm_order);

                    % Only use complete bins
                    usable_samples = n_bins * n_samples;
                    data = data(1:usable_samples, :, :, :);

                    % Reshape to group samples
                    new_sz = [n_samples, n_bins, sz(perm_order(2:end))];
                    reshaped = reshape(data, new_sz);

                    % Calculate mean along first dimension (now our target dimension)
                    valid_counts = sum(~isnan(reshaped), 1);
                    totals = sum(reshaped, 1, 'omitnan');
                    mean_data = squeeze(totals ./ valid_counts);

                    % Reshape back to original dimensionality
                    final_sz = sz;
                    final_sz(dim_idx) = n_bins;

                    % Permute back to original dimension order
                    inv_perm_order(perm_order) = 1:length(perm_order);
                    ds_out.(current_field).data = permute(reshape(mean_data, [n_bins sz(perm_order(2:end))]), inv_perm_order);

                    % Update coordinates for this dimension if they exist
                    if isfield(ds.(current_field), 'coords')
                        coord_fields = fieldnames(ds.(current_field).coords);
                        for k = 1:length(coord_fields)
                            coord_field = coord_fields{k};
                            if contains(lower(coord_field), lower(dim_to_find))
                                coord_data = ds.(current_field).coords.(coord_field);
                                ds_out.(current_field).coords.(coord_field) = coord_data(1:n_samples:usable_samples);
                            end
                        end
                    end
                end
            end
        end
    end

    % Update the main coordinate system
    if isfield(ds, 'coords')
        coord_fields = fieldnames(ds.coords);
        for i = 1:length(coord_fields)
            if contains(lower(coord_fields{i}), lower(dim_to_find))
                coord_data = ds.coords.(coord_fields{i});
                usable_samples = floor(length(coord_data)/n_samples) * n_samples;
                ds_out.coords.(coord_fields{i}) = coord_data(1:n_samples:usable_samples);
            end
        end
    end
end

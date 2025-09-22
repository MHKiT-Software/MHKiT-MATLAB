function reshaped_data = reshape_for_ensemble(data, n_bin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reshape dolfyn time series data into ensemble bins
%
% This function takes continuous time series data and reshapes it into
% ensemble bins, enabling calculation of statistics within each ensemble.
%
% Parameters
% ------------
%   data: array
%       Input data with time as the last dimension
%       Shape: (..., time)
%   n_bin: integer  
%       Number of time points per ensemble bin
%
% Returns
% ---------
%   reshaped_data: array
%       Reshaped data with ensemble structure
%       Shape: (..., n_ensembles, n_bin)
%       where n_ensembles = floor(time_length / n_bin)
%
% Notes
% -----
% - Incomplete ensembles at the end are discarded
% - Time dimension must be last dimension of input data
% - Equivalent to MHKiT-Python: data[..., :length*n_bin].reshape(..., length, n_bin)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Validate inputs
    if ~isnumeric(data)
        error('MHKiT:dolfyn:reshape_for_ensemble', 'Data must be numeric array');
    end
    
    if ~isscalar(n_bin) || n_bin <= 0 || mod(n_bin, 1) ~= 0
        error('MHKiT:dolfyn:reshape_for_ensemble', 'n_bin must be positive integer');
    end
    
    % Get dimensions
    sz = size(data);
    time_length = sz(end);
    
    % Calculate number of complete ensembles
    n_ensembles = floor(time_length / n_bin);
    
    if n_ensembles == 0
        error('MHKiT:dolfyn:reshape_for_ensemble', ...
            'n_bin (%d) is larger than time dimension (%d)', n_bin, time_length);
    end
    
    % Calculate usable time points (discard incomplete ensemble)
    usable_length = n_ensembles * n_bin;
    
    % Extract only complete ensembles
    if length(sz) == 1
        % 1D case: [time] -> [n_ensembles, n_bin]
        data_truncated = data(1:usable_length);
        reshaped_data = reshape(data_truncated, [n_bin, n_ensembles])';
    elseif length(sz) == 2
        % 2D case: [dim1, time] -> [dim1, n_ensembles, n_bin] 
        data_truncated = data(:, 1:usable_length);
        reshaped_data = reshape(data_truncated, [sz(1), n_bin, n_ensembles]);
        reshaped_data = permute(reshaped_data, [1, 3, 2]); % [dim1, n_ensembles, n_bin]
    else
        % Multi-dimensional case: [..., time] -> [..., n_ensembles, n_bin]
        % Use linear indexing for the last dimension
        data_truncated = data;
        sz_new = sz;
        sz_new(end) = usable_length;
        
        % Reshape to collapse all but last dimension
        data_2d = reshape(data_truncated, [prod(sz(1:end-1)), sz(end)]);
        data_2d_truncated = data_2d(:, 1:usable_length);
        
        % Reshape to ensemble format
        data_ensemble = reshape(data_2d_truncated, [prod(sz(1:end-1)), n_bin, n_ensembles]);
        data_ensemble = permute(data_ensemble, [1, 3, 2]); % [prod(dims), n_ensembles, n_bin]
        
        % Reshape back to original dimensions plus ensemble structure
        final_shape = [sz(1:end-1), n_ensembles, n_bin];
        reshaped_data = reshape(data_ensemble, final_shape);
    end
    
    % Warn if data was discarded
    if usable_length < time_length
        warning('MHKiT:dolfyn:reshape_for_ensemble', ...
            'Discarded %d time points to create %d complete ensembles of %d points each', ...
            time_length - usable_length, n_ensembles, n_bin);
    end

end

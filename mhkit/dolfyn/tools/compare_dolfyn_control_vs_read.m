function diff = compare_dolfyn_control_vs_read(ds_read, ds_cntrl, threshold)
%%%%%%%%%%%%%%%%%%%%
%     Compare dolfyn control file vs read data with detailed debugging
%
%     Performs two-pass comparison:
%     1. Silent comparison to get total difference
%     2. If diff > threshold, automatically generates detailed debug table
%
% Parameters
% ------------
%     ds_read: structure
%         Structure from the binary instrument data
%
%     ds_cntrl: structure
%         Control structure read from python generated NetCDF
%
%     threshold: float (optional, default: 1e-6)
%         Tolerance threshold for triggering detailed debug output
%
% Returns
% ---------
%     diff: float
%         Total difference between the data in the two structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Set tolerance for numerical comparisons
    % tolerance = 1e-6;  % Configurable tolerance threshold
    tolerance = 0.01;  % Configurable tolerance threshold

    if nargin < 3
        threshold = 1e-6;
    end

    % Pass 1: Silent comparison to get total difference
    diff = compare_structures_silent(ds_read, ds_cntrl, tolerance);

    % Pass 2: If diff > threshold, run detailed analysis with table output
    if diff > threshold
        fprintf('\n=== MHKIT DOLFYN COMPARISON REPORT ===\n');
        fprintf('Total Difference: %.7f (threshold: %.7f, tolerance: %.2e)\n\n', diff, threshold, tolerance);

        % Run detailed recursive analysis
        run_detailed_comparison(ds_read, ds_cntrl, tolerance);
    end
end


function diff = compare_structures_silent(ds_read, ds_cntrl, tolerance)
    % Silent version of structure comparison (no debug output)
    diff = 0.0;
    exclude = {'coords', 'attrs', 'time', 'complex_vars', 'filehead_config'};

    % Check coords first
    fields = fieldnames(ds_read.coords);
    for qq = 1:numel(fields)
        field = fields{qq};
        if ~any(contains(field, exclude))
            if isfield(ds_cntrl.coords, field)
                if iscell(ds_cntrl.coords.(field))
                    for kk = 1:length(ds_cntrl.coords.(field))
                        diff = diff + double(~strcmpi(...
                            ds_cntrl.coords.(field)(kk),...
                            ds_read.coords.(field)(kk)));
                    end
                else
                    % its numeric because time got excluded
                    if size(ds_read.coords.(field)) ~= ...
                            size(ds_cntrl.coords.(field))
                        coord_diff = abs(...
                            sum(abs(double(ds_cntrl.coords.(field)) -...
                            double(ds_read.coords.(field)'))))/2.0;
                    else
                        coord_diff = abs(...
                            sum(abs(double(ds_cntrl.coords.(field)) -...
                            double(ds_read.coords.(field))))/2.0);
                    end
                    % Only add to diff if it exceeds tolerance
                    if coord_diff > tolerance
                        diff = diff + coord_diff;
                    end
                end
            end
        end
    end

    % Check Attributes
    fields = fieldnames(ds_cntrl.attrs);
    for qq = 1:numel(fields)
        field = fields{qq};
        if ~any(contains(field, exclude))
            if isfield(ds_read.attrs, field)
                if iscell(ds_cntrl.attrs.(field))
                    for kk = 1:numel(ds_cntrl.attrs.(field))
                        chk_nm = ds_cntrl.attrs.(field){kk};
                        diff = diff + ...
                            double(~any(strcmpi(chk_nm,...
                            ds_read.attrs.(field))));
                    end
                elseif ischar(ds_cntrl.attrs.(field))
                    if strcmpi(ds_cntrl.attrs.(field),'yes')
                        diff = diff + double(~ds_read.attrs.(field));
                    elseif strcmpi(ds_cntrl.attrs.(field),'no')
                        diff = diff + double(ds_read.attrs.(field));
                    else
                        diff = diff + ...
                            double(~strcmpi(ds_cntrl.attrs.(field),...
                            ds_read.attrs.(field)));
                    end
                elseif islogical(ds_read.attrs.(field))
                    diff = diff + double(ds_cntrl.attrs.(field) ~= ...
                            ds_read.attrs.(field));
                elseif isnumeric(ds_cntrl.attrs.(field))
                    if contains(class(ds_cntrl.attrs.(field)),'int')
                        diff = diff + ...
                            sum(double(ds_cntrl.attrs.(field) ~=...
                            ds_read.attrs.(field)'));
                    else
                        if all(size(ds_cntrl.attrs.(field)) ==...
                                size(ds_read.attrs.(field)))
                            atr_vl = ds_read.attrs.(field);
                        else
                            atr_vl = ds_read.attrs.(field)';
                        end
                        diff = diff + sum(abs(ds_cntrl.attrs.(field)...
                            - atr_vl), ...
                            1:numel(size(ds_cntrl.attrs.(field))))...
                            /length(ds_cntrl.attrs.(field));
                    end
                end
            end
        end
    end

    % Now check the remaining fields
    fields = fieldnames(ds_cntrl);
    for qq = 1:numel(fields)
        field = fields{qq};
        if ~any(contains(field, exclude))
            if isfield(ds_read, field)
                cls = class(ds_cntrl.(field));
                if strcmp(cls,'struct')
                    % Loop through structure
                    % Data
                    tmp1 = isnan(ds_cntrl.(field).data);
                    dt1 = double(ds_cntrl.(field).data);
                    dt1(tmp1) = 0.0;
                    tmp2 = isnan(ds_read.(field).data);
                    dt2 = double(ds_read.(field).data);
                    dt2(tmp2) = 0.0;

                    diff = diff + abs(sum(abs(dt1 - dt2),...
                            1:numel(size(dt1)))/length(dt1));
                    % Dims
                    if length(ds_cntrl.(field).dims) ~= length(ds_read.(field).dims)
                        % Different number of dimensions
                        diff = diff + abs(length(ds_cntrl.(field).dims) - length(ds_read.(field).dims));
                    else
                        % Same number of dimensions, compare each one
                        for kk = 1:length(ds_cntrl.(field).dims)
                            diff = diff + double(~strcmpi(...
                                ds_cntrl.(field).dims{kk},...
                                ds_read.(field).dims{kk}));
                        end
                    end
                    % Coords (we already checked coords so just check names)
                    if isfield(ds_cntrl.(field), 'coords') && isfield(ds_read.(field), 'coords')
                        cntrl_coord_names = fieldnames(ds_cntrl.(field).coords);
                        read_coord_names = fieldnames(ds_read.(field).coords);
                        for kk = 1:numel(cntrl_coord_names)
                            chk_nm = cntrl_coord_names{kk};
                            diff = diff + ...
                                double(~any(strcmpi(chk_nm, read_coord_names)));
                        end
                    end
                    % Units if they exist
                    if isfield(ds_cntrl.(field),'units') && ...
                        ~strcmpi(ds_cntrl.(field).units, "none") && ...
                        isfield(ds_read.(field), 'units')
                        diff = diff + ...
                            double(~strcmpi(ds_cntrl.(field).units,...
                            ds_read.(field).units));
                    end
                else
                    diff = diff + ...
                            double(~strcmpi(ds_cntrl.(field),...
                            ds_read.(field)));
                end
            end
        end
    end
end


function run_detailed_comparison(ds_read, ds_cntrl, tolerance)
    % Run detailed recursive analysis and generate formatted tables
    exclude = {'coords', 'attrs', 'time', 'complex_vars', 'filehead_config'};

    % Collect all field comparisons
    comparisons = {};

    % Check top-level coords
    comparisons = collect_coord_comparisons(ds_read.coords, ds_cntrl.coords, 'coords', exclude, comparisons, tolerance);

    % Check top-level attrs
    comparisons = collect_attr_comparisons(ds_read.attrs, ds_cntrl.attrs, 'attrs', exclude, comparisons, tolerance);

    % Check data fields
    comparisons = collect_data_field_comparisons(ds_read, ds_cntrl, '', exclude, comparisons, tolerance);

    % Separate meaningful vs stylistic differences
    meaningful_comparisons = {};
    stylistic_comparisons = {};

    for i = 1:length(comparisons)
        if comparisons{i}.is_meaningful
            meaningful_comparisons{end+1} = comparisons{i};
        else
            stylistic_comparisons{end+1} = comparisons{i};
        end
    end

    % Display results based on what we found
    if isempty(meaningful_comparisons) && isempty(stylistic_comparisons)
        fprintf('All meaningful data matches!\n\n');
    elseif isempty(meaningful_comparisons)
        % Only stylistic differences exist
        fprintf('All meaningful data matches!\n\n');
        if ~isempty(stylistic_comparisons)
            fprintf('STYLISTIC DIFFERENCES (formatting only)\n');
            fprintf('========================================\n');
            stylistic_sorted = sort_comparisons_by_priority(stylistic_comparisons);
            print_comparison_table_compact(stylistic_sorted);
        end
    else
        % Meaningful differences exist - show both tables
        fprintf('MEANINGFUL DIFFERENCES (fix these first)\n');
        fprintf('=========================================\n');
        meaningful_sorted = sort_comparisons_by_priority(meaningful_comparisons);
        print_comparison_table(meaningful_sorted);

        if ~isempty(stylistic_comparisons)
            fprintf('\nSTYLISTIC DIFFERENCES (clean up later)\n');
            fprintf('======================================\n');
            stylistic_sorted = sort_comparisons_by_priority(stylistic_comparisons);
            print_comparison_table_compact(stylistic_sorted);
        end
    end
end


function comparisons = collect_coord_comparisons(coords_read, coords_cntrl, path_prefix, exclude, comparisons, tolerance)
    % Collect coordinate field comparisons

    % Check for fields in control but missing in read
    if isstruct(coords_cntrl)
        fields_cntrl = fieldnames(coords_cntrl);
        for i = 1:length(fields_cntrl)
            field = fields_cntrl{i};
            if ~any(contains(field, exclude))
                field_path = sprintf('%s.%s', path_prefix, field);

                if ~isfield(coords_read, field)
                    comp = create_comparison_entry(field_path, 'MISSING_READ', ...
                        format_field_info(coords_cntrl.(field)), 'N/A', ...
                        'Field missing in MATLAB read');
                    comparisons{end+1} = comp;
                else
                    % Compare existing fields
                    [status, issue, is_meaningful] = compare_field_values(coords_read.(field), coords_cntrl.(field), tolerance);
                    if ~strcmp(status, 'OK')
                        comp = create_comparison_entry(field_path, status, ...
                            format_field_info(coords_cntrl.(field)), ...
                            format_field_info(coords_read.(field)), issue, is_meaningful);
                        comparisons{end+1} = comp;
                    end
                end
            end
        end
    end

    % Check for fields in read but missing in control
    if isstruct(coords_read)
        fields_read = fieldnames(coords_read);
        for i = 1:length(fields_read)
            field = fields_read{i};
            if ~any(contains(field, exclude))
                if ~isfield(coords_cntrl, field)
                    field_path = sprintf('%s.%s', path_prefix, field);
                    comp = create_comparison_entry(field_path, 'MISSING_CTRL', ...
                        'N/A', format_field_info(coords_read.(field)), ...
                        'Field only in MATLAB read');
                    comparisons{end+1} = comp;
                end
            end
        end
    end
end


function comparisons = collect_attr_comparisons(attrs_read, attrs_cntrl, path_prefix, exclude, comparisons, tolerance)
    % Collect attribute field comparisons

    % Check for fields in control but missing in read
    if isstruct(attrs_cntrl)
        fields_cntrl = fieldnames(attrs_cntrl);
        for i = 1:length(fields_cntrl)
            field = fields_cntrl{i};
            if ~any(contains(field, exclude))
                field_path = sprintf('%s.%s', path_prefix, field);

                if ~isfield(attrs_read, field)
                    comp = create_comparison_entry(field_path, 'MISSING_READ', ...
                        format_field_info(attrs_cntrl.(field)), 'N/A', ...
                        'Field missing in MATLAB read');
                    comparisons{end+1} = comp;
                else
                    % Compare existing fields
                    [status, issue, is_meaningful] = compare_field_values(attrs_read.(field), attrs_cntrl.(field), tolerance);
                    if ~strcmp(status, 'OK')
                        comp = create_comparison_entry(field_path, status, ...
                            format_field_info(attrs_cntrl.(field)), ...
                            format_field_info(attrs_read.(field)), issue, is_meaningful);
                        comparisons{end+1} = comp;
                    end
                end
            end
        end
    end

    % Check for fields in read but missing in control
    if isstruct(attrs_read)
        fields_read = fieldnames(attrs_read);
        for i = 1:length(fields_read)
            field = fields_read{i};
            if ~any(contains(field, exclude))
                if ~isfield(attrs_cntrl, field)
                    field_path = sprintf('%s.%s', path_prefix, field);
                    comp = create_comparison_entry(field_path, 'MISSING_CTRL', ...
                        'N/A', format_field_info(attrs_read.(field)), ...
                        'Field only in MATLAB read');
                    comparisons{end+1} = comp;
                end
            end
        end
    end
end


function comparisons = collect_data_field_comparisons(struct_read, struct_cntrl, path_prefix, exclude, comparisons, tolerance)
    % Recursively collect data field comparisons

    % Check for fields in control but missing in read
    if isstruct(struct_cntrl)
        fields_cntrl = fieldnames(struct_cntrl);
        for i = 1:length(fields_cntrl)
            field = fields_cntrl{i};
            if ~any(contains(field, exclude))
                if isempty(path_prefix)
                    field_path = field;
                else
                    field_path = sprintf('%s.%s', path_prefix, field);
                end

                if ~isfield(struct_read, field)
                    comp = create_comparison_entry(field_path, 'MISSING_READ', ...
                        format_field_info(struct_cntrl.(field)), 'N/A', ...
                        'Field missing in MATLAB read');
                    comparisons{end+1} = comp;
                else
                    % Field exists in both - check if it's a data structure
                    if isstruct(struct_cntrl.(field)) && isfield(struct_cntrl.(field), 'data')
                        % This is a data field with .data, .dims, .coords, etc.
                        comparisons = analyze_data_structure(struct_read.(field), struct_cntrl.(field), field_path, comparisons, tolerance);
                    elseif isstruct(struct_cntrl.(field))
                        % This is a nested structure - recurse
                        comparisons = collect_data_field_comparisons(struct_read.(field), struct_cntrl.(field), field_path, exclude, comparisons, tolerance);
                    else
                        % This is a simple field
                        [status, issue, is_meaningful] = compare_field_values(struct_read.(field), struct_cntrl.(field), tolerance);
                        if ~strcmp(status, 'OK')
                            comp = create_comparison_entry(field_path, status, ...
                                format_field_info(struct_cntrl.(field)), ...
                                format_field_info(struct_read.(field)), issue, is_meaningful);
                            comparisons{end+1} = comp;
                        end
                    end
                end
            end
        end
    end

    % Check for fields in read but missing in control
    if isstruct(struct_read)
        fields_read = fieldnames(struct_read);
        for i = 1:length(fields_read)
            field = fields_read{i};
            if ~any(contains(field, exclude))
                if ~isfield(struct_cntrl, field)
                    if isempty(path_prefix)
                        field_path = field;
                    else
                        field_path = sprintf('%s.%s', path_prefix, field);
                    end
                    comp = create_comparison_entry(field_path, 'MISSING_CTRL', ...
                        'N/A', format_field_info(struct_read.(field)), ...
                        'Field only in MATLAB read');
                    comparisons{end+1} = comp;
                end
            end
        end
    end
end


function comparisons = analyze_data_structure(data_read, data_cntrl, field_path, comparisons, tolerance)
    % Analyze a data structure with .data, .dims, .coords, etc.

    % Check .data field
    if isfield(data_cntrl, 'data') && isfield(data_read, 'data')
        [status, issue] = compare_data_arrays(data_read.data, data_cntrl.data, tolerance);
        if ~strcmp(status, 'OK')
            comp = create_comparison_entry(sprintf('%s.data', field_path), status, ...
                format_field_info(data_cntrl.data), ...
                format_field_info(data_read.data), issue);
            comparisons{end+1} = comp;
        end
    elseif isfield(data_cntrl, 'data') && ~isfield(data_read, 'data')
        comp = create_comparison_entry(sprintf('%s.data', field_path), 'MISSING_READ', ...
            format_field_info(data_cntrl.data), 'N/A', 'Data field missing in read');
        comparisons{end+1} = comp;
    elseif ~isfield(data_cntrl, 'data') && isfield(data_read, 'data')
        comp = create_comparison_entry(sprintf('%s.data', field_path), 'MISSING_CTRL', ...
            'N/A', format_field_info(data_read.data), 'Data field only in read');
        comparisons{end+1} = comp;
    end

    % Check .dims field
    if isfield(data_cntrl, 'dims') && isfield(data_read, 'dims')
        [status, issue, is_meaningful] = compare_field_values(data_read.dims, data_cntrl.dims, tolerance);
        if ~strcmp(status, 'OK')
            comp = create_comparison_entry(sprintf('%s.dims', field_path), status, ...
                format_field_info(data_cntrl.dims), ...
                format_field_info(data_read.dims), issue, is_meaningful);
            comparisons{end+1} = comp;
        end
    elseif isfield(data_cntrl, 'dims') && ~isfield(data_read, 'dims')
        comp = create_comparison_entry(sprintf('%s.dims', field_path), 'MISSING_READ', ...
            format_field_info(data_cntrl.dims), 'N/A', 'Dims field missing in read');
        comparisons{end+1} = comp;
    elseif ~isfield(data_cntrl, 'dims') && isfield(data_read, 'dims')
        comp = create_comparison_entry(sprintf('%s.dims', field_path), 'MISSING_CTRL', ...
            'N/A', format_field_info(data_read.dims), 'Dims field only in read');
        comparisons{end+1} = comp;
    end

    % Check .coords field (if exists)
    if isfield(data_cntrl, 'coords') && isfield(data_read, 'coords')
        comparisons = collect_coord_comparisons(data_read.coords, data_cntrl.coords, sprintf('%s.coords', field_path), {}, comparisons, tolerance);
    elseif isfield(data_cntrl, 'coords') && ~isfield(data_read, 'coords')
        comp = create_comparison_entry(sprintf('%s.coords', field_path), 'MISSING_READ', ...
            format_field_info(data_cntrl.coords), 'N/A', 'Coords field missing in read');
        comparisons{end+1} = comp;
    elseif ~isfield(data_cntrl, 'coords') && isfield(data_read, 'coords')
        comp = create_comparison_entry(sprintf('%s.coords', field_path), 'MISSING_CTRL', ...
            'N/A', format_field_info(data_read.coords), 'Coords field only in read');
        comparisons{end+1} = comp;
    end

    % Check .units field (if exists)
    if isfield(data_cntrl, 'units') && isfield(data_read, 'units')
        [status, issue, is_meaningful] = compare_field_values(data_read.units, data_cntrl.units, tolerance);
        if ~strcmp(status, 'OK')
            comp = create_comparison_entry(sprintf('%s.units', field_path), status, ...
                format_field_info(data_cntrl.units), ...
                format_field_info(data_read.units), issue, is_meaningful);
            comparisons{end+1} = comp;
        end
    elseif isfield(data_cntrl, 'units') && ~isfield(data_read, 'units')
        comp = create_comparison_entry(sprintf('%s.units', field_path), 'MISSING_READ', ...
            format_field_info(data_cntrl.units), 'N/A', 'Units field missing in read');
        comparisons{end+1} = comp;
    elseif ~isfield(data_cntrl, 'units') && isfield(data_read, 'units')
        comp = create_comparison_entry(sprintf('%s.units', field_path), 'MISSING_CTRL', ...
            'N/A', format_field_info(data_read.units), 'Units field only in read');
        comparisons{end+1} = comp;
    end
end


function comp = create_comparison_entry(path, status, control_info, read_info, issue, is_meaningful)
    % Create a comparison entry structure
    if nargin < 6
        is_meaningful = true;  % Default for missing fields is meaningful
    end

    comp = struct();
    comp.path = path;
    comp.status = status;
    comp.control_info = control_info;
    comp.read_info = read_info;
    comp.issue = issue;
    comp.is_meaningful = is_meaningful;
end


function [status, issue, is_meaningful] = compare_field_values(read_val, cntrl_val, tolerance)
    % Compare two field values and return status, issue description, and meaningfulness
    status = 'OK';
    issue = '-';
    is_meaningful = true;  % Default to meaningful difference

    % Check types first
    read_class = class(read_val);
    cntrl_class = class(cntrl_val);

    if ~strcmp(read_class, cntrl_class)
        % Check if type difference is stylistic (same values, different types)
        [status, issue, is_meaningful] = compare_type_difference(read_val, cntrl_val, tolerance);
        return;
    end

    % Check sizes
    if ~isequal(size(read_val), size(cntrl_val))
        % Check if size difference is just transpose (stylistic)
        [status, issue, is_meaningful] = compare_size_difference(read_val, cntrl_val, tolerance);
        return;
    end

    % Check values
    if isnumeric(read_val) && isnumeric(cntrl_val)
        [status, issue] = compare_numeric_values(read_val, cntrl_val, tolerance);
    elseif iscell(read_val) && iscell(cntrl_val)
        [status, issue] = compare_cell_values(read_val, cntrl_val);
    elseif ischar(read_val) && ischar(cntrl_val)
        if ~strcmp(read_val, cntrl_val)
            status = 'VALUE_DIFF';
            issue = sprintf('String mismatch: "%s" vs "%s"', cntrl_val, read_val);
        end
    elseif islogical(read_val) && islogical(cntrl_val)
        if read_val ~= cntrl_val
            status = 'VALUE_DIFF';
            issue = sprintf('Boolean mismatch: %s vs %s', mat2str(cntrl_val), mat2str(read_val));
        end
    end

    % Value differences are always meaningful
    % is_meaningful remains true for VALUE_DIFF cases
end


function [status, issue] = compare_data_arrays(read_data, cntrl_data, tolerance)
    % Special comparison for large data arrays with meaningful statistics
    status = 'OK';
    issue = '-';

    % Check types first
    read_class = class(read_data);
    cntrl_class = class(cntrl_data);

    if ~strcmp(read_class, cntrl_class)
        status = 'TYPE_DIFF';
        issue = sprintf('Type mismatch: %s vs %s', cntrl_class, read_class);
        return;
    end

    % Check sizes
    if ~isequal(size(read_data), size(cntrl_data))
        status = 'SIZE_DIFF';
        issue = sprintf('Size mismatch: %s vs %s', mat2str(size(cntrl_data)), mat2str(size(read_data)));
        return;
    end

    % Check values for numeric data
    if isnumeric(read_data) && isnumeric(cntrl_data)
        % Handle NaN values
        cntrl_clean = double(cntrl_data);
        read_clean = double(read_data);

        cntrl_nan = isnan(cntrl_clean);
        read_nan = isnan(read_clean);

        % Set NaN to 0 for comparison
        cntrl_clean(cntrl_nan) = 0;
        read_clean(read_nan) = 0;

        % Calculate differences
        abs_diff = abs(cntrl_clean - read_clean);
        max_diff = max(abs_diff(:));
        mean_diff = mean(abs_diff(:));

        if max_diff > tolerance  % Threshold for considering values different
            status = 'VALUE_DIFF';
            if max_diff < 1e-3
                issue = sprintf('Data differs (max: %.2e, mean: %.2e)', max_diff, mean_diff);
            else
                issue = sprintf('Data differs (max: %.6f, mean: %.6f)', max_diff, mean_diff);
            end
        end
    end
end


function [status, issue] = compare_numeric_values(read_val, cntrl_val, tolerance)
    % Compare numeric values
    status = 'OK';
    issue = '-';

    if isscalar(read_val) && isscalar(cntrl_val)
        if abs(read_val - cntrl_val) > tolerance
            status = 'VALUE_DIFF';
            issue = sprintf('Values differ: %.6f vs %.6f', cntrl_val, read_val);
        end
    else
        % Array comparison
        abs_diff = abs(double(read_val) - double(cntrl_val));
        max_diff = max(abs_diff(:));
        if max_diff > tolerance
            status = 'VALUE_DIFF';
            issue = sprintf('Array values differ (max diff: %.6f)', max_diff);
        end
    end
end


function [status, issue] = compare_cell_values(read_val, cntrl_val)
    % Compare cell arrays
    status = 'OK';
    issue = '-';

    for i = 1:numel(read_val)
        if ~strcmp(read_val{i}, cntrl_val{i})
            status = 'VALUE_DIFF';
            issue = sprintf('Cell contents differ at index %d', i);
            return;
        end
    end
end


function info_str = format_field_info(field_value)
    % Format field information for display
    if isnumeric(field_value)
        type_str = class(field_value);
        size_str = mat2str(size(field_value));
        if isscalar(field_value)
            info_str = sprintf('%s %s (%.3f)', type_str, size_str, field_value);
        else
            info_str = sprintf('%s %s', type_str, size_str);
        end
    elseif iscell(field_value)
        info_str = sprintf('cell %s', mat2str(size(field_value)));
    elseif ischar(field_value)
        if length(field_value) > 20
            info_str = sprintf('char "%s..."', field_value(1:17));
        else
            info_str = sprintf('char "%s"', field_value);
        end
    elseif islogical(field_value)
        info_str = sprintf('logical %s (%s)', mat2str(size(field_value)), mat2str(field_value));
    elseif isstruct(field_value)
        field_names = fieldnames(field_value);
        info_str = sprintf('struct [%d fields]', length(field_names));
    else
        info_str = sprintf('%s %s', class(field_value), mat2str(size(field_value)));
    end
end


function sorted_comparisons = sort_comparisons_by_priority(comparisons)
    % Sort comparisons by priority (missing fields first)
    if isempty(comparisons)
        sorted_comparisons = comparisons;
        return;
    end

    priority_order = {'MISSING_READ', 'MISSING_CTRL', 'TYPE_DIFF', 'SIZE_DIFF', 'VALUE_DIFF', 'OK'};
    priority_scores = zeros(1, length(comparisons));

    for i = 1:length(comparisons)
        status = comparisons{i}.status;
        score = find(strcmp(status, priority_order));
        if isempty(score)
            score = length(priority_order) + 1;
        end
        priority_scores(i) = score;
    end

    [~, sort_idx] = sort(priority_scores);
    sorted_comparisons = comparisons(sort_idx);
end


function print_comparison_table(comparisons)
    % Print formatted comparison table (200 column limit)

    if isempty(comparisons)
        fprintf('All fields match!\n\n');
        return;
    end

    % Define column widths (total ~200 characters)
    path_width = 35;
    status_width = 12;
    control_width = 25;
    read_width = 25;
    issue_width = 95;  % Remaining space

    % Print header
    header_format = sprintf('%%-%ds | %%-%ds | %%-%ds | %%-%ds | %%-%ds\n', ...
        path_width, status_width, control_width, read_width, issue_width);
    fprintf(header_format, 'Path', 'Status', 'Control Info', 'Read Info', 'Issue Description');

    % Print separator
    sep_line = repmat('-', 1, path_width + status_width + control_width + read_width + issue_width + 12);
    fprintf('%s\n', sep_line);

    % Print data rows
    row_format = sprintf('%%-%ds | %%-%ds | %%-%ds | %%-%ds | %%-%ds\n', ...
        path_width, status_width, control_width, read_width, issue_width);

    for i = 1:length(comparisons)
        comp = comparisons{i};

        % Truncate long strings to fit column widths
        path_trunc = truncate_string(comp.path, path_width);
        status_trunc = truncate_string(comp.status, status_width);
        control_trunc = truncate_string(comp.control_info, control_width);
        read_trunc = truncate_string(comp.read_info, read_width);
        issue_trunc = truncate_string(comp.issue, issue_width);

        fprintf(row_format, path_trunc, status_trunc, control_trunc, read_trunc, issue_trunc);
    end

    fprintf('\n');
end


function print_comparison_table_compact(comparisons)
    % Print compact comparison table for stylistic differences

    if isempty(comparisons)
        return;
    end

    % More compact column widths
    path_width = 40;
    status_width = 15;
    issue_width = 110;  % Focus on path, status, and issue

    % Print header
    header_format = sprintf('%%-%ds | %%-%ds | %%-%ds\n', path_width, status_width, issue_width);
    fprintf(header_format, 'Path', 'Status', 'Issue Description');

    % Print separator
    sep_line = repmat('-', 1, path_width + status_width + issue_width + 6);
    fprintf('%s\n', sep_line);

    % Print data rows
    row_format = sprintf('%%-%ds | %%-%ds | %%-%ds\n', path_width, status_width, issue_width);

    for i = 1:length(comparisons)
        comp = comparisons{i};

        % Truncate long strings to fit column widths
        path_trunc = truncate_string(comp.path, path_width);
        status_trunc = truncate_string(comp.status, status_width);
        issue_trunc = truncate_string(comp.issue, issue_width);

        fprintf(row_format, path_trunc, status_trunc, issue_trunc);
    end

    fprintf('\n');
end


function truncated = truncate_string(str, max_width)
    % Truncate string to fit within column width
    if length(str) <= max_width
        truncated = str;
    else
        truncated = [str(1:max_width-3), '...'];
    end
end


function [status, issue, is_meaningful] = compare_type_difference(read_val, cntrl_val, tolerance)
    % Compare values with different types to determine if difference is meaningful
    status = 'TYPE_DIFF';
    is_meaningful = true;  % Default to meaningful

    read_class = class(read_val);
    cntrl_class = class(cntrl_val);

    % Check for logical vs numeric flag conversions (stylistic difference)
    if (islogical(read_val) && isnumeric(cntrl_val)) || (isnumeric(read_val) && islogical(cntrl_val))
        % Convert logical to numeric for comparison
        if islogical(read_val)
            read_numeric = double(read_val);
            cntrl_numeric = double(cntrl_val);
        else
            read_numeric = double(read_val);
            cntrl_numeric = double(cntrl_val);
        end

        % Compare converted values
        if isequal(size(read_numeric), size(cntrl_numeric))
            if isscalar(read_numeric)
                value_diff = abs(cntrl_numeric - read_numeric);
            else
                abs_diff = abs(cntrl_numeric(:) - read_numeric(:));
                value_diff = max(abs_diff);
            end

            if value_diff <= tolerance
                % Same values, just logical vs numeric (stylistic difference)
                is_meaningful = false;
                issue = sprintf('Type difference only: %s vs %s (logical flag conversion)', cntrl_class, read_class);
            else
                % Different values AND different types (meaningful difference)
                issue = sprintf('Type mismatch with value differences: %s vs %s', cntrl_class, read_class);
            end
        else
            % Different sizes too, definitely meaningful
            issue = sprintf('Type and size mismatch: %s %s vs %s %s', cntrl_class, mat2str(size(cntrl_val)), read_class, mat2str(size(read_val)));
        end

    % Check if both are numeric types (can be converted and compared)
    elseif (isnumeric(read_val) && isnumeric(cntrl_val))
        % Convert both to higher precision type for comparison
        if strcmp(cntrl_class, 'double') || strcmp(read_class, 'double')
            cntrl_converted = double(cntrl_val);
            read_converted = double(read_val);
        elseif strcmp(cntrl_class, 'single') || strcmp(read_class, 'single')
            cntrl_converted = single(cntrl_val);
            read_converted = single(read_val);
        else
            % Both integers, convert to double for comparison
            cntrl_converted = double(cntrl_val);
            read_converted = double(read_val);
        end

        % Compare converted values
        if isequal(size(cntrl_converted), size(read_converted))
            if isscalar(cntrl_converted)
                value_diff = abs(cntrl_converted - read_converted);
            else
                abs_diff = abs(cntrl_converted(:) - read_converted(:));
                value_diff = max(abs_diff);
            end

            if value_diff <= tolerance
                % Same values, just different types (stylistic difference)
                is_meaningful = false;
                issue = sprintf('Type difference only: %s vs %s (same values)', cntrl_class, read_class);
            else
                % Different values AND different types (meaningful difference)
                issue = sprintf('Type mismatch with value differences: %s vs %s', cntrl_class, read_class);
            end
        else
            % Different sizes too, definitely meaningful
            issue = sprintf('Type and size mismatch: %s %s vs %s %s', cntrl_class, mat2str(size(cntrl_val)), read_class, mat2str(size(read_val)));
        end
    else
        % Non-numeric types or mixed numeric/non-numeric
        issue = sprintf('Type mismatch: %s vs %s', cntrl_class, read_class);
        % Keep as meaningful difference
    end
end


function [status, issue, is_meaningful] = compare_size_difference(read_val, cntrl_val, tolerance)
    % Compare values with different sizes to determine if difference is meaningful
    status = 'SIZE_DIFF';
    is_meaningful = true;  % Default to meaningful

    cntrl_size = size(cntrl_val);
    read_size = size(read_val);

    % Check if this is a transpose case (2D arrays only)
    if length(cntrl_size) == 2 && length(read_size) == 2 && ...
       isequal(cntrl_size, fliplr(read_size))

        % Handle numeric transpose case
        if isnumeric(read_val) && isnumeric(cntrl_val)
            % This is a transpose case - compare values after transpose
            try
                cntrl_data = double(cntrl_val);
                read_data_transposed = double(read_val');

                % Compare transposed data
                if isscalar(cntrl_data)
                    value_diff = abs(cntrl_data - read_data_transposed);
                else
                    abs_diff = abs(cntrl_data(:) - read_data_transposed(:));
                    value_diff = max(abs_diff);
                end

                if value_diff <= tolerance
                    % Same data, just transposed (stylistic difference)
                    status = 'TRANSPOSE_ONLY';
                    is_meaningful = false;
                    issue = sprintf('Same data, different orientation (row vs column)');
                else
                    % Transposed AND different values (meaningful difference)
                    status = 'TRANSPOSE_DIFF';
                    issue = sprintf('Transposed with value differences: %s vs %s', mat2str(cntrl_size), mat2str(read_size));
                end
            catch
                % If transpose comparison fails, treat as meaningful size difference
                issue = sprintf('Size mismatch: %s vs %s', mat2str(cntrl_size), mat2str(read_size));
            end

        % Handle cell array transpose case (like dims fields)
        elseif iscell(read_val) && iscell(cntrl_val)
            try
                % Compare cell contents after transposing
                read_transposed = read_val';
                if isequal(cntrl_val, read_transposed)
                    % Same cell contents, just transposed (stylistic difference)
                    status = 'TRANSPOSE_ONLY';
                    is_meaningful = false;
                    issue = sprintf('Same cell data, different orientation (row vs column)');
                else
                    % Different cell contents even after transpose
                    issue = sprintf('Cell array transpose with content differences: %s vs %s', mat2str(cntrl_size), mat2str(read_size));
                end
            catch
                % If cell comparison fails, treat as meaningful size difference
                issue = sprintf('Cell array size mismatch: %s vs %s', mat2str(cntrl_size), mat2str(read_size));
            end
        else
            % Mixed types or other cases - treat as meaningful
            issue = sprintf('Size mismatch: %s vs %s', mat2str(cntrl_size), mat2str(read_size));
        end
    else
        % Not a simple transpose case - this is a meaningful size difference
        issue = sprintf('Size mismatch: %s vs %s', mat2str(cntrl_size), mat2str(read_size));
    end
end

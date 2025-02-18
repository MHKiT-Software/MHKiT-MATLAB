function output = typecast_from_mhkit_python(mhkit_python_output)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Convert Python output to appropriate MATLAB type
%
% Parameters
% ------------
%     mhkit_python_output: Python structure
%         Python output structure containing type information and data
%
% Returns
% ---------
%     output: structure
%         Converted data in MATLAB-compatible format containing:
%         - data: The converted data in appropriate MATLAB format
%         - index: (Optional) Structure containing index information:
%             * name: Name of the index
%             * data: The index values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % First convert Python data to native types using Python utility
    py_output = py.mhkit_python_utils.type_conversion.convert_to_matlab_compatible(mhkit_python_output);

    % Get type, data, and index from converted output
    dataType = char(py_output{'type'});
    data = py_output{'data'};
    index = py_output{'index'};

    % Initialize output structure
    output = struct();

    % Handle data conversion based on type
    switch dataType
        case '<class ''numpy.ndarray''>'
            output.data = double(data);

        case 'int'
            output.data = int64(data);

        case 'double'
            output.data = double(data);

        case 'scalar'
            if isa(data, 'py.dict')
                temp = struct();
                keys = cell(py.list(data.keys()));
                for i = 1:length(keys)
                    key = char(keys{i});
                    value = data{keys{i}};
                    if isscalar(value)
                        if mod(value, 1) == 0
                            temp.(key) = int64(value);
                        else
                            temp.(key) = double(value);
                        end
                    else
                        % Handle numpy array
                        temp.(key) = double(value);
                    end
                end
                output.data = temp;
            else
                output.data = double(data);
            end

        case {'scalar_int', 'array_int'}
            if isa(data, 'py.dict')
                temp = struct();
                keys = cell(py.list(data.keys()));
                for i = 1:length(keys)
                    key = char(keys{i});
                    value = data{keys{i}};
                    % Direct conversion of numpy arrays
                    temp.(key) = int64(value);
                end
                output.data = temp;
            else
                output.data = int64(data);
            end

        case {'scalar_double', 'array_double'}
            keys = cell(py.list(data.keys()));
            output.data = struct();
            output.columns = cellfun(@char, keys, 'UniformOutput', false);

            if length(keys) == 1
                output.data = double(data{keys{1}});
            else
                for i = 1:length(keys)
                    key = char(keys{i});
                    sanitized_key = sanitize_fieldname(key);
                    % Direct conversion of numpy arrays
                    output.data.(sanitized_key) = double(data{keys{i}});
                end
            end
            % if isa(data, 'py.dict')
            % else
            %     output.data = double(data);
            % end

        case 'mixed'
            if ~isa(data, 'py.dict')
                error('Mixed type data should be a dictionary');
            end
            temp = struct();
            keys = cell(py.list(data.keys()));
            for i = 1:length(keys)
                key = char(keys{i});
                value = data{keys{i}};
                if isscalar(value)
                    if mod(value, 1) == 0
                        temp.(key) = int64(value);
                    else
                        temp.(key) = double(value);
                    end
                else
                    % Direct conversion of numpy arrays
                    temp.(key) = double(value);
                end
            end
            output.data = temp;

        case 'capture_length_matrix'
            output = struct();
            % Direct conversion of numpy arrays
            output.x_centers = double(data{"x_centers"});
            output.y_centers = double(data{"y_centers"});
            output.capture_length_matrix = double(data{"capture_length_matrix"});

            return

        otherwise
            error('Unsupported type: %s', dataType);
    end

    output.index = struct();

    % Add index information if available
    if ~isempty(index) && ~isprop(index, 'empty') && ~ismethod(index, 'empty')
        if isa(index, 'py.dict')
            output.index.name = char(index{'name'});
            output.index.data = double(index{'data'});
        end
    end
end

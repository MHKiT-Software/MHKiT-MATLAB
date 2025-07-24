function write_netcdf(ds, filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Write Dolfyn data set to NetCDF data structure.
%     Updated for MATLAB R2022b and later versions.
%
% Parameters
% ------------
%     ds: structure
%         Structure from the binary instrument data
%
%     filename: string
%         Filename of NetCDF file to read.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make sure that ds is a structure
    if ~isstruct(ds)
        ME = MException('MATLAB:write_netcdf','ds must be a structure');
        throw(ME);
    end

    % Make sure that ds contains the dolfyn fields
    if ~isfield(ds,'coords') || ~isfield(ds,'attrs') || ...
            ~isfield(ds,'time')
        ME = MException('MATLAB:write_netcdf',['The provided data ' ...
            'structure does not appear to have been created by dolfyn']);
        throw(ME);
    end

    % Ensure proper file extension
    if ~endsWith(filename, ".nc")
        [filepath, name, ~] = fileparts(filename);
        filename = fullfile(filepath, name + ".nc");
    end

    % Check for file existence and handle overwrite
    overwriteFile = true;
    if exist(filename, 'file')
        promptMessage = sprintf('This file already exists:\n%s\nDo you want to overwrite it?', filename);
        titleBarCaption = 'Overwrite?';
        buttonText = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
        if strcmpi(buttonText, 'No')
            overwriteFile = false;
        end
    end

    if overwriteFile
        if exist(filename, 'file')
            delete(filename);
        end

        % Process coordinates
        fields = string(fieldnames(ds.coords));
        dim_map = containers.Map(fields, 1:numel(fields));

        for idx = 1:numel(fields)
            key = fields(idx);
            n = numel(ds.coords.(key));

            if iscell(ds.coords.(key))
                data = string(ds.coords.(key));
                type = 'string';
            else
                data = ds.coords.(key);
                type = class(data);
            end

            % Create dimension and variable
            nccreate(filename, key, 'Dimensions', {key, n}, ...
                     'Datatype', type, 'Format', 'netcdf4');
            ncwrite(filename, key, data);
        end

        % Process variables
        exclude = ["coords", "coord_sys", "attrs", "time", "hdwtime_gps"];
        fields = string(fieldnames(ds));

        for idx = 1:numel(fields)
            key = fields(idx);
            if ~any(exclude == key)
                dim_fields = string(fieldnames(ds.(key).coords));
                dimensions = cell(1, numel(dim_fields)*2);

                for k = 1:numel(dim_fields)
                    dimensions{(k-1)*2 + 1} = dim_fields(k);
                    dimensions{(k-1)*2 + 2} = numel(ds.(key).coords.(dim_fields(k)));
                end

                data = ds.(key).data;
                if islogical(data)
                    data = int32(data);
                end

                % Create and write variable
                nccreate(filename, key, 'Dimensions', dimensions, ...
                         'Datatype', class(data), 'FillValue', nan, ...
                         'Format', 'netcdf4');
                ncwrite(filename, key, squeeze(data));

                % Handle units
                if isfield(ds.(key), 'units')
                    units = ds.(key).units;
                else
                    units = 'None';
                end
                ncwriteatt(filename, key, 'units', units);
            end
        end

        % Process attributes
        fields = string(fieldnames(ds.attrs));
        for idx = 1:numel(fields)
            key = fields(idx);
            if ~isstruct(ds.attrs.(key))
                attr_val = ds.attrs.(key);
                if iscell(attr_val)
                    attr_val = string(attr_val);
                elseif islogical(attr_val)
                    attr_val = int32(attr_val);
                end
                ncwriteatt(filename, '/', key, attr_val);
            end
        end

        % Add empty complex_vars attribute for Python compatibility
        ncwriteatt(filename, '/', 'complex_vars', []);
    end
end

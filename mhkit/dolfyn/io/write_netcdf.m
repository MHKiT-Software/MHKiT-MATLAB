function write_netcdf(ds, filename)
%%%%%%%%%%%%%%%%%%%%
%     Write Dolfyn data set to NetCDF data structure.
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

    % Make sure that ds contains the dolyn fields
    if ~isfield(ds,'coords') || ~isfield(ds,'attrs') || ...
            ~isfield(ds,'time')
        ME = MException('MATLAB:write_netcdf',['The provided data ' ...
            'structure does not appear to have been created by dolfyn']);
        throw(ME);
    end

    overwriteFile = true; % Default to overwriting or creating new file.
    if exist(filename, 'file')
        % Ask user if they want to overwrite the file.
        promptMessage = sprintf(['This file already exists:\n%s\nDo you ' ...
          'want to overwrite it?'], filename);
        titleBarCaption = 'Overwrite?';
        buttonText = questdlg(promptMessage, ...
          titleBarCaption, 'Yes', 'No', 'Yes');
        if strcmpi(buttonText, 'No')
        % User does not want to overwrite. 
        % Set flag to not do the write.
        overwriteFile = false;
        end
    end
    if overwriteFile
        % File does not exist yet, or the user wants to overwrite
        % an existing file.
        if exist(filename, 'file')
            delete(filename);
        end
        % Loop through coords and create netcdf dimensions
        fields = fieldnames(ds.coords);
        for qq = 1:numel(fields)
            key = fields{qq};
            n = numel(ds.coords.(key));            
            if iscell(ds.coords.(key))
                temp = convertCharsToStrings(ds.coords.(key));
                type = 'string';
            else
                type = class(ds.coords.(key));
            end
            nccreate(filename, key,'Dimensions',{key, n},...
                    'Datatype', type,...
                    'Format','netcdf4');
            if iscell(ds.coords.(key))
                ncwrite(filename, key, temp);
            elseif numel(ds.coords.(key)) == 4
                ncwrite(filename, key, ds.coords.(key)(1:3));
            else
                ncwrite(filename, key, ds.coords.(key));
            end
        end        
        
        % List of fields to exclude from the variable write
        exclude = {'coords', 'coord_sys', 'attrs', 'time', 'hdwtime_gps'};
        % Loop through the ds fields and create the variable in the netcdf
        % file then write the data.
        fields = fieldnames(ds);
        for qq = 1:numel(fields)
            key = fields{qq};
            if strcmp(key,'orientmat')
                debug = 1;
            end
            if ~any(strcmp(exclude,key))
                dim_fields = fieldnames(ds.(key).coords);
                dimensions = cell(1,numel(dim_fields)*2);
                for kk = 1:numel(dim_fields)
                    dimensions{(kk-1)*2 + 1} = dim_fields{kk};
                    n = numel(ds.(key).coords.(dim_fields{kk}));                    
                    dimensions{(kk-1)*2 + 2} = n;                        
                end
                % dimensions cell can now be used to create the netcdf
                % variable
                nccreate(filename, key,'Dimensions',dimensions,...
                    'Datatype', class(ds.(key).data),...
                    'FillValue', nan, 'Format','netcdf4');
                % Now that the variable exists we can write the data to it 
                out_data = squeeze(ds.(key).data);
                ncwrite(filename, key, out_data);
                % add the units
                if isfield(ds.(key), 'units')
                    ncwriteatt(filename, key, "units", [ds.(key).units]);
                else
                    ncwriteatt(filename, key, "units", 'None');
                end
            end
        end 
   
        % Attributes
        fields = fieldnames(ds.attrs);
        for qq = 1:numel(fields)
            key = fields{qq};
            if ~isstruct(ds.attrs.(key))
                if iscell(ds.attrs.(key))
                    ncwriteatt(filename, '/', key, ...
                        convertCharsToStrings(ds.attrs.(key)));
                elseif islogical(ds.attrs.(key))
                    ncwriteatt(filename, '/', key, int32(ds.attrs.(key)));
                else
                    ncwriteatt(filename, '/', key, ds.attrs.(key));
                end
            end
        end
        % Not currently functional but needed for python counterpart
        ncwriteatt(filename, '/', 'complex_vars', []);
    end

end


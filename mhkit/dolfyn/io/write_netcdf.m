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
    if ~isMATLABReleaseOlderThan("R2021b")
        if ~endsWith(filename, ".h5")
            temp = extractBetween(filename,1,strfind(filename, "."));
            filename = join([temp{1},"h5"],"");
        end
    end
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

        if ~isMATLABReleaseOlderThan("R2021b")
            % Older versions of Matlab did not allow strings to be
            % written to netcdf so we have to use h5
            fun_map = struct(  ...
                'create',{@create_h5},...
                'write',{@write_h5},...
                'attribute',{@attr_h5}); 
        else
            fun_map = struct(  ...
                'create',{@create_nc},...
                'write',{@write_nc},...
                'attribute',{@attr_nc});
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
            feval(fun_map.create, filename, key, {key, n}, type, true);
            if iscell(ds.coords.(key))
                feval(fun_map.write, filename, key, temp, true);            
            else
                feval(fun_map.write, filename, key, ds.coords.(key), true);                
            end
        end        
        
        % List of fields to exclude from the variable write
        exclude = {'coords', 'coord_sys', 'attrs', 'time', 'hdwtime_gps'};
        % Loop through the ds fields and create the variable in the netcdf
        % file then write the data.
        fields = fieldnames(ds);
        for qq = 1:numel(fields)
            key = fields{qq};            
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
                if islogical(ds.(key).data)
                    ds.(key).data = int32(ds.(key).data);
                end
                feval(fun_map.create, filename, key, dimensions,...
                    class(ds.(key).data), false);
                % Now that the variable exists we can write the data to it 
                out_data = squeeze(ds.(key).data);
                feval(fun_map.write, filename, key, out_data, false);
                % add the units
                if isfield(ds.(key), 'units')
                    feval(fun_map.attribute, filename, key, ...
                        [ds.(key).units], false);
                else
                    feval(fun_map.attribute, filename, key, 'None', false);
                end
            end
        end 
   
        % Attributes
        fields = fieldnames(ds.attrs);
        for qq = 1:numel(fields)
            key = fields{qq};
            if ~isstruct(ds.attrs.(key))
                if iscell(ds.attrs.(key))
                    feval(fun_map.attribute, filename, key, ...
                        convertCharsToStrings(ds.attrs.(key)), true);
                elseif islogical(ds.attrs.(key))
                    feval(fun_map.attribute, filename, key, ...
                        int32(ds.attrs.(key)), true);
                else
                    feval(fun_map.attribute, filename, key, ...
                        ds.attrs.(key), true);
                end
            end
        end
        % Not currently functional but needed for python counterpart
        feval(fun_map.attribute, filename, 'complex_vars', [], true);            
    end   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NetCDF Functions 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function create_nc(filename, key, dim_cell, type, coords)
        if coords
            nccreate(filename, key,'Dimensions',dim_cell,...
			    'Datatype', type,'Format','netcdf4');
        else
            nccreate(filename, key,'Dimensions',dim_cell,...
			'Datatype', type,'FillValue', nan, 'Format','netcdf4');
        end
    end

    function write_nc(filename, key, data, ~)
        ncwrite(filename, key, data)
    end

    function attr_nc(filename, key, attribute, attrs)
        if attrs
            ncwriteatt(filename,'/', key, attribute);
        else
            ncwriteatt(filename, key, "units", attribute);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % H5 Functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function create_h5(filename, key, dim_cell, type, coords)
        if coords
            loc = join(['/coords/',key],'');
            size = dim_cell{2};
        else
            loc = join(['/',key],'');
            size = zeros(1,numel(dim_cell)/2);
            dims = '';
            for i = 1:numel(dim_cell)/2
                size(i) = dim_cell{i*2};
                dims = join([dims,"/",dim_cell{(i-1)*2 + 1}],'');
            end
        end
        h5create(filename, loc, size,'Datatype', type);
        if ~coords
            h5writeatt(filename, join(['/',key]), "dims", dims);
        end
    end

    function write_h5(filename, key, data, coords)
        if coords
            h5write(filename, join(['/coords/',key],''), data);
        else
            h5write(filename, join(['/',key],''), data);
        end
    end

    function attr_h5(filename, key, attribute, attrs)
        if attrs
            h5writeatt(filename,'/', key, attribute);
        else
            h5writeatt(filename, join(['/',key]), "units", attribute);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


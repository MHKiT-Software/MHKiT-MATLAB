function ds = read_h5(filename)
%%%%%%%%%%%%%%%%%%%%
%     Read H5 data structure.
%     
% Parameters
% ------------
%     filename: string
%         Filename of NetCDF file to read.
%
% Returns
% ---------
%     ds: structure 
%         Structure from the binary instrument data
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % check to see if the filename input is a string
    if ~ischar(filename)
        ME = MException('MATLAB:read_h5',['filename must be a ' ...
            'character string']);
        throw(ME);
    end
    
    % check to see if the file exists
    if ~isfile(filename)
        ME = MException('MATLAB:read_h5','file does not exist');
        throw(ME);
    end

    ds = struct();

    % Get NetCDF info to populate the variable names    
    info = h5info(filename); 

    % Check if the file was written by netcdf 
    if strcmp(info.Attributes(1).Name,'_NCProperties')
        nc_file = true;
    else
        nc_file = false;
    end
    coord_vals = {};
    coord_keys = [];
  
    % Loop through the variables once to get the coords
    var_size = numel(info.Datasets);
    for qq = 1:var_size
        is_coord = false;
        name = info.Datasets(qq).Name;        
        path = join(['/',name],"");
        dtype = info.Datasets(qq).Datatype.Class;
        attrs = info.Datasets(qq).Attributes;
        for kk=1:numel(attrs)
            if strcmpi(attrs(kk).Name,'CLASS') && ...
                    strcmpi(attrs(kk).Value,'dimension_scale')            
                % variable is a coordinate
                is_coord = true;
            elseif strcmpi(attrs(kk).Name,'_Netcdf4Dimid')
                coord_vals(end+1) = {name};
                coord_keys(end+1) = attrs(kk).Value;
            end
        end
        if is_coord            
            if strcmp(dtype, 'H5T_STRING')
                % if its a string then we read a string array
                ds.coords.(name) = ...
                    convertStringsToChars(h5read(filename,path));
            else                
                if strcmpi(name, 'x*')
                    coord_vals(end) = {'x_star'};
                    ds.coords.x_star = h5read(filename,path);
                else
                    ds.coords.(name) = h5read(filename,path);
                end                
            end 
        end
    end

    dim_map = containers.Map(coord_keys,coord_vals);

    % Loop through the variables again to get the remaining data
    for qq = 1:var_size
        name = info.Datasets(qq).Name;
        path = join(['/',name],"");
        attrs = info.Datasets(qq).Attributes;
        sz = info.Datasets(qq).ChunkSize;
        units = 'None';
        if ~any(strcmp(name,coord_vals))
            % variable goes in the main structure field
            for i= 1:numel(attrs)
                if strcmpi(attrs(i).Name,'_Netcdf4Coordinates')
                    dimensions = cell(numel(attrs(i).Value),1);
                    if nc_file
                        count = 1;
                        for jj = numel(attrs(i).Value):-1:1                       
                            dimensions(count) = ...
                                {dim_map(attrs(i).Value(jj))};
                            count = count + 1;
                        end 
                    else
                        for jj = 1:numel(attrs(i).Value)                        
                            dimensions(jj) = {dim_map(attrs(i).Value(jj))};
                        end     
                    end
                elseif strcmpi(attrs(i).Name,'units') 
                    if iscell(attrs(i).Value)
                        units = attrs(i).Value{1};
                    else
                        units = attrs(i).Value;
                    end
                end
            end
            if numel(sz) == 1
                % no modifications needed (the read function does it)
                ds.(name).data = h5read(filename,path);                               
            elseif numel(sz) == 2
                if contains(name, 'orientmat') || ...
                        contains(name,'inst2head_rotmat')
                    % no modifications needed 
                    ds.(name).data = h5read(filename,path);
                else
                    % Need to reshape the data
                    temp_dat = h5read(filename,path); 
                    tmp_shape = size(temp_dat);
                    tmp_shape = [tmp_shape(1),1,tmp_shape(2)];
                    temp_dat = reshape(temp_dat,tmp_shape);
                    ds.(name).data = temp_dat;
                end             
            else
                % Need to reshape the data
                temp_dat = h5read(filename,path); 
                tmp_shape = size(temp_dat);
                tmp_shape = [tmp_shape(1),1,tmp_shape(2:3)];
                temp_dat = reshape(temp_dat,tmp_shape);
                ds.(name).data = temp_dat;                                
            end
            ds.(name).dims = dimensions;
            for kk = 1:numel(dimensions)
                ds.(name).coords.(dimensions{kk}) = ...
                    ds.coords.(dimensions{kk});
            end 
            ds.(name).units = units;
        end
    end

    % Finally grab the attributes
    for qq = 1:numel(info.Attributes)
        name = info.Attributes(qq).Name;
        if startsWith(name,'_')
            continue
        end
        value = info.Attributes(qq).Value;
        ds.attrs.(name) = value;
    end

    % Corrections 
    if isfield(ds.attrs,'rotate_vars')
        ds.attrs.rotate_vars = cellstr(ds.attrs.rotate_vars);
    end
    if isfield(ds, 'orientation_down')
        ds.orientation_down.data = logical(ds.orientation_down.data);
    end

    ds.coord_sys = ds.attrs.coord_sys;
    ds.time = ds.coords.time;

end
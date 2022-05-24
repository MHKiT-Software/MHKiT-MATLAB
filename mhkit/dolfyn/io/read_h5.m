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
    
    % Start with collecting the coords  
    dim_size = numel(info.Groups.Datasets);
%     coords = cell(dim_size,1);
    for qq = 1:dim_size
        name = info.Groups.Datasets(qq).Name;
        path = join(['/coords/', name], ''); 
        dtype = info.Groups.Datasets(qq).Datatype.Class;
        if strcmp(dtype, 'H5T_STRING')
            % if its a string then we read a string array
            ds.coords.(name) = ...
                convertStringsToChars(h5read(filename,path));
        else                
            if strcmpi(name, 'x*')
                ds.coords.x_star = h5read(filename,path);
            else
                ds.coords.(name) = h5read(filename,path);
            end                
        end  
    end
    
    % Loop through the data variables 
    var_size = numel(info.Datasets);
    for qq = 1:var_size
        name = info.Datasets(qq).Name;
        path = join(['/',name],"");
        sz = info.Datasets(qq).Dataspace.Size;
        for i= 1:numel(info.Datasets(qq).Attributes)
            if strcmpi(info.Datasets(qq).Attributes(i).Name,'dims')
                dimensions = info.Datasets(qq).Attributes(i).Value;
                dimensions = split(dimensions,'/');
                dimensions(1) = [];
            elseif strcmpi(info.Datasets(qq).Attributes(i).Name,'units')        
                attrs = info.Datasets(qq).Attributes(i).Value;
            end
        end
        
        % variable goes in the main structure field
        if numel(sz) == 1
            % no modifications needed (the read function does it)
            ds.(name).data = h5read(filename,path);
            ds.(name).dims = cell(numel(dimensions),1);
            for kk = 1:numel(dimensions)
                ds.(name).dims{kk} = dimensions{kk};
                ds.(name).coords.(dimensions{kk}) = ...
                    ds.coords.(dimensions{kk});
            end
            if ~isempty(attrs)
                if ischar(attrs)
                    ds.(name).units = attrs;
                else
                    ds.(name).units = attrs{1};
                end
            end
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
            ds.(name).dims = cell(numel(dimensions),1);
            for kk = 1:numel(dimensions)
                if strcmpi(dimensions{kk}, 'x*')
                    ds.(name).dims{kk} = 'x_star';
                    ds.(name).coords.x_star = ds.coords.x_star;
                else
                    ds.(name).dims{kk} = dimensions{kk};
                    ds.(name).coords.(dimensions{kk}) = ...
                        ds.coords.(dimensions{kk});
                end
            end
            if ~isempty(attrs)
                if ischar(attrs)
                    ds.(name).units = attrs;
                else
                    ds.(name).units = attrs{1};
                end
            end
        else
            % Need to reshape the data
            temp_dat = h5read(filename,path); 
            tmp_shape = size(temp_dat);
            tmp_shape = [tmp_shape(1),1,tmp_shape(2:3)];
            temp_dat = reshape(temp_dat,tmp_shape);
            ds.(name).data = temp_dat;
            ds.(name).dims = cell(numel(dimensions),1);
            for kk = 1:numel(dimensions)
                if strcmpi(dimensions{kk}, 'x*')
                    ds.(name).dims{kk} = 'x_star';
                else
                    ds.(name).dims{kk} = dimensions{kk};
                end
                ds.(name).coords.(dimensions{kk}) = ...
                    ds.coords.(dimensions{kk});
            end
            if ~isempty(attrs)
                if ischar(attrs)
                    ds.(name).units = attrs;
                else
                    ds.(name).units = attrs{1};
                end
            end               
        end
        
    end

    % Finally grab the attributes
    for qq = 1:numel(info.Attributes)
        name = info.Attributes(qq).Name;
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
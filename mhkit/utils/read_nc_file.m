function ds = read_nc_file_var(filename,varargin)
%%%%%%%%%%%%%%%%%%%%
%     Read NetCDF data structure.
%     
% Parameters
% ------------
%     filename: string
%         Filename of NetCDF file to read.
%
%     vnms: cell array of characters (optional)
%         variable names to read.
%         if the variable is in a group, make sure to include group name,
%         i.e., {'GROUP_NAME1/VAR_NAME1','GROUP_NAME2/VAR_NAME2',...}
%         read all variables under '/' if vnms not given.
%
% Returns
% ---------
%     ds: structure 
%         Structure from the binary instrument data
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % check to see if the filename input is a string
    if ~ischar(filename) && ~isstring(filename)
        ME = MException('MATLAB:read_nc_file_var',['filename must be a ' ...
            'character string']);
        throw(ME);
    % check to see if the file exists
    elseif ~isfile(filename)
        ME = MException('MATLAB:read_netcdf','file does not exist');
        throw(ME);
    % check MATLAB version & if it is h5 file
    elseif isMATLABReleaseOlderThan("R2021b") || ...
            endsWith(filename, ".h5")
        ds = read_h5(filename); 
        return
    end

    finfo = ncinfo(filename);
    % assign variable names to vnms
    if nargin > 1
        % use user specified input
        vnms = varargin{1};
    elseif ~isempty(finfo.Variables)
        % read all variables under '/' 
        vnms = {finfo.Variables.Name};
    else
        % return ERROR if 
        % no variable specified and no Variables under '/'
        ME = MException('MATLAB:read_nc_file',['no variable available' ...
            ' to read']);
        throw(ME);
    end

    % initiate output struct()
    ds = struct();
    % loop through variable names (vnms) to read data and attributes
    for ivar=1:numel(vnms)
        name = vnms{ivar};tmplst = split(name,'/');
        vname = tmplst{end};
        if nargin > 1
            % use user specified input
            vinfo = ncinfo(filename,name);
        else
            % read all variables under '/' 
            vinfo = finfo.Variables(ivar);
        end
        % check if vname is valid to be a field name of struct()
        % if not, convert it to a valid field name
        vname = check_name(vname);
        % read in metadata and assign it to the variable
        ds.(vname).data = ncread(filename,name);
        % assign dims name to the variable
        if ~isempty(vinfo.Dimensions)
            ds.(vname).dims = {vinfo.Dimensions.Name};
        else
            ds.(vname).dims = {};
        end
        % assign FillValue to the variable
        if ~isnumeric(vinfo.FillValue)
            ds.(vname).FillValue = str2double(vinfo.FillValue);
        else
            ds.(vname).FillValue = vinfo.FillValue;
        end
        ds.(vname).attrs = struct();
        % if no arributes skip
        if isempty(vinfo.Attributes)
            continue
        end
        % otherwise, loop through to 
        % assign Attributes & check _FillValue
        attrnames = {vinfo.Attributes.Name};
        for iattr = 1:numel(attrnames)
            aname = attrnames{iattr};
            % _FillValue in data attributes:
            if strcmp(aname,'_FillValue') 
                % assign a new FillValue only if it was NaN
                if isnan(ds.(vname).FillValue)
                    ds.(vname).FillValue = ...
                        vinfo.Attributes(iattr).Value;
                end
                continue
            else
                aname = check_name(aname);
                ds.(vname).attrs.(aname) = ...
                    vinfo.Attributes(iattr).Value;
            end
        end
    end
    % assign global attributes of the file
    if ~isempty(finfo.Attributes)
        attrnames = {finfo.Attributes.Name};
        for iattr = 1:numel(attrnames)
            aname = attrnames{iattr};
            aname = check_name(aname);
            ds.attrs.(aname) = finfo.Attributes(iattr).Value;
        end
    end
end

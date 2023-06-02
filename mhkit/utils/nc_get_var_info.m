function res = nc_get_var_info(filename,vinfo,varpath,res)
%%%%%%%%%%%%%%%%%%%%
%     Append NetCDF variable data and info to an existring structure (res).
%     
% Parameters
% ------------
%   filename: string
%       Filename of NetCDF file to read.
%   vinfo: 
%       NetCDF variable schema, a structure.
%   varpath:
%       Variable path in the NetCDF file, 
%       E.G., 'GROUP_NAME/SUBGROUP_NAME/VAR_NAME'
%   res: 
%       Structure to hold the variable info
%       
%
% Returns
% ---------
%   res: structure that holds the variable info
%       fields will include variables from nc file 
%       e.g., res.(varname1), res.(varname2), ...
%       res.(varname) includes fields: 
%           res.(varname).Name: full name of the variable
%           res.(varname).Data: metadata 
%           res.(varname).Dims: dimension names
%           res.(varname).FillValue: V.FillValue or 
%               V.Attributes{'_FillValue'} if
%               the former is not given.
%           res.(varname).Attrs: V.Attributes except for 
%               V.Attributes{'_FillValue'}
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    vname = vinfo.Name;
    % check & convert vname
    vname = check_name(vname);
    res.(vname).Name = vinfo.Name;
    
    % assign data to the variable
    res.(vname).Data = ncread(filename, varpath);

    % assign dims to the variable
    if ~isempty(vinfo.Dimensions)
        res.(vname).Dims = {vinfo.Dimensions.Name};
    else
        res.(vname).Dims = {};
    end

    % assign FillValue to the variable
    if ~isnumeric(vinfo.FillValue)
        res.(vname).FillValue = ...
            str2double(vinfo.FillValue);
    else
        res.(vname).FillValue = vinfo.FillValue;
    end
    
    res.(vname).Attrs = struct();
    % assign attrs to the variable & check _FillValue
    if ~isempty(vinfo.Attributes)
        attrnames = {vinfo.Attributes.Name};
        for iattr = 1:numel(attrnames)
            aname = attrnames{iattr};
            % _FillValue in data attributes:
            if strcmp(aname,'_FillValue') 
                % assign a new FillValue only if it was NaN
                if isnan(res.(vname).FillValue)
                    res.(vname).FillValue = ...
                        vinfo.Attributes(iattr).Value;
                end
                continue
            else
                aname = check_name(aname);
                res.(vname).Attrs.(aname) = ...
                    vinfo.Attributes(iattr).Value;
            end
        end
    end
end
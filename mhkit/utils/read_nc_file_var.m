function ds = read_nc_file_var(filename,vnms,out_opt)
%%%%%%%%%%%%%%%%%%%%
%     Read NetCDF data structure.
%     
% Parameters
% ------------
%     filename: string
%         Filename of NetCDF file to read.
%
%     vnms: cell array of characters
%         variable names to read.
%         if the variable is in a group, make sure to include group name,
%         i.e., {'GROUP_NAME1/VAR_NAME1','GROUP_NAME2/VAR_NAME2',...}
%
%     out_opt: 0 if output as structure, ds.(vname) & ds.attrs.
%              1 if output as structure array, ds(N).
%
% Returns
% ---------
%     ds: structure or structure array
%       1. structure that has fields:
%          1.1 variable names from nc file 
%       e.g., res.(varname1), res.(varname2), ...
%       res.(varname) includes fields: 
%           res.(varname).Name: full name of the variable
%           res.(varname).FillValue: V.FillValue or 
%               V.Attributes{'_FillValue'} if
%               the former is not given.
%           res.(varname).Data: metadata 
%           res.(varname).Dims: dimension names
%           res.(varname).Attrs: V.Attributes except for 
%               V.Attributes{'_FillValue'}
%           1.2 global attributes from nc file
%       
%       2. structure array, and each structure has the same fields as
%       res.(varname) stated above.
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % precheck:
    nc_file_precheck(filename);

    finfo = ncinfo(filename);
    % if want a struct array as output
    if out_opt ==1
        % Use data struct array to hold:
        N = numel(vnms);
        ds(N)=struct('Name',[],'Data',[], ...
            'Dims',[],'FillValue',[],'Attrs',[]);
        % loop through variable names (vnms) to read data and attributes
        for ivar=1:numel(vnms)
            name = vnms{ivar};
            vinfo = ncinfo(filename,name);
            ds(ivar).Name = vinfo.Name;
            ds(ivar).FillValue = str2double(vinfo.FillValue);
            ds(ivar).Data = ncread(filename, name);
            ds(ivar).Dims = {vinfo.Dimensions.Name};
            ds(ivar).Attrs = vinfo.Attributes;
        end
        return
    end
    % if want a struct
    ds = struct();
    % loop through variable names (vnms) to read data and attributes
    for ivar=1:numel(vnms)
        name = vnms{ivar};
        vinfo = ncinfo(filename,name);
        ds = nc_get_var_info(filename,vinfo,name,ds);
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

function ds = read_nc_file_var(filename,vnms,opt_out)
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
%     opt_out: 0 will output as structure, ds.(vname) & ds.attrs.
%              1 will output as structure array, ds(N).
%   
%   Note 1. Variables in root group ('/') can be extracted 
%    directly by 'VARNAME'.
%   Note 2. To find the group info of a variable, other than reading 
%    the provided user manual of the corresponding NetCDF data product, 
%    one can use ncdisp(filename) in MATLAB which shows variable info 
%    for all groups and subgroups. 
%   An example ncdisp(filename.nc) output can be:
%   ...
%       /GrpName/SubGrpName/
%       Variables: 
%           ...
%           VAR_I_Need
%               Size: ...
%               Dimensions: ...
%               ...
%   To get VAR_I_Need from filename.nc, one can run this function as:
% 
%    ds = read_nc_file_var('filename.nc',...
%           {'/GrpName/SubGrpName/VAR_I_Need',...
%           '/GrpName/SubGrpName/Another_VAR_I_Need'},...
%           opt_out=0);
%   to get a struct.
% 
%    OR: 
%    ds = read_nc_file_var('filename.nc',...
%           {'/GrpName/SubGrpName/VAR_I_Need',...
%           '/GrpName/SubGrpName/Another_VAR_I_Need'},...
%           opt_out=1);
%   to get a struct array.
%   
%   Note 3. If needs to read all variables within a Group, one can use ncinfo():
%       
%    ginfo = ncinfo('filename.nc','/GrpName/SubGrpName/');
%    vnms = {ginfo.Variables.Name};
%    vnms = strcat('GrpName/SubGrpName/',vnms);
%    ds = read_nc_file_var('filename.nc',vnms,opt_out=0);%struct 
%    %ds = read_nc_file_var('filename.nc',vnms,opt_out=1);%struct array
%
% Returns
% ------------
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
%       2. struct array, and each struct has the same fields as
%       res.(varname) stated above.
%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % precheck:
    nc_file_precheck(filename);

    finfo = ncinfo(filename);
    % if want a struct array as output
    if opt_out ==1
        % Use data struct array to hold:
        N = numel(vnms);
        ds(N)=struct('Name',[],'Data',[], ...
            'Dims',[],'FillValue',[],'Attrs',[]);
        % loop through variable names (vnms) to read data and attributes
        for ivar=1:numel(vnms)
            name = vnms{ivar};
            vinfo = ncinfo(filename,name);
            ds(ivar).Name = vinfo.Name;
            if isstring(vinfo.FillValue)
                ds(ivar).FillValue = str2double(vinfo.FillValue);
            else
                ds(ivar).FillValue = vinfo.FillValue;
            end
            ds(ivar).Data = ncread(filename, name);
            ds(ivar).Attrs = vinfo.Attributes;
            if ~isempty(vinfo.Dimensions)
                ds(ivar).Dims = {vinfo.Dimensions.Name};
            end
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

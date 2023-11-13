function data = read_nc_file(fname,varargin)

%%%%%%%%%%%%%%%%%%%%
%     Read NetCDF data into a MATLAB struct().
%
% Parameters
% ------------
%   filename: string
%       Filename of NetCDF file to read.
%   vnms: cell array of characters (optional)
%       Variable names to read.
%       Read all variables if vnms not given. Otherwise, read specified
%       variables.
%       If the variable is in a group, make sure to include group name,
%       i.e., {'GROUP_NAME1/VAR_NAME1',
%              'GROUP_NAME2/SUBGROUP_NAME2/VAR_NAME2',...}
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
%    ds = read_nc_file('filename.nc',...
%           {'/GrpName/SubGrpName/VAR_I_Need'});
%
%   Note 3. If needs to read all variables within a Group, one can use ncinfo():
%
%    ginfo = ncinfo('filename.nc','/GrpName/SubGrpName/');
%    vnms = {ginfo.Variables.Name};
%    vnms = strcat('GrpName/SubGrpName/',vnms);
%    ds = read_nc_file('filename.nc',vnms);
%
% Returns
% ------------
%     ds: struct
%   1. If vnms not specified, ds has fields:
%       groups: struct that contains names of child groups as fields
%       LongName: path to current group
%       Attributes: Attributes of current group
%       Variables: struct that contains each variable as a struct and
%           a list of all variables within this group.
%           e.g., Variables.(varname1), Variables.(varname2), ...
%         Variables.(varname) includes fields:
%           Name: full name of the variable
%           Data: metadata
%           Dims: dimension names
%           FillValue: V.FillValue or V.Attributes{'_FillValue'} if
%               the former is not given.
%           Attrs: V.Attributes except for V.Attributes{'_FillValue'}
%   2. If vnms specified, ds has fields: ds.(varname).
%       Each ds.(varname) contains fields same as Variables.(varname)
%       described above.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % precheck:
    nc_file_precheck(fname);

    % if variables specified:
    if nargin > 1
        vnms = varargin{1};
        data = read_nc_file_var(fname,vnms,0);
        return
    end

    % if variables not specified:
    data = struct();
    finfo = ncinfo(fname);
    data = nc_read_file_group(fname,finfo,data,'');
    data = data.root;

    function ds = nc_read_file_group(filename,ginfo,ds,long_name)
        gname = ginfo.Name;
        if strcmp(gname,'/')
            gname = 'root';
            ds.(gname) = struct();
            ds.(gname).LongName = '';
        else
            gname = check_name(gname);
            ds.(gname) = struct();
            ds.(gname).LongName = strcat(long_name,'/',ginfo.Name);
        end
        if ~isempty(ginfo.Variables)
            ds.(gname).Variables = struct();
            ds.(gname).Variables.AllVarNames = {ginfo.Variables.Name};
            for ivar = 1:numel(ginfo.Variables)
                vinfo = ginfo.Variables(ivar);
                varpath = strcat(ds.(gname).LongName,'/',vinfo.Name);
                ds.(gname).Variables = nc_get_var_info(filename,vinfo,...
               varpath,ds.(gname).Variables);
            end

        end

        if ~isempty(ginfo.Attributes)
            ds.(gname).Attributes = ginfo.Attributes;
        end

        if ~isempty(ginfo.Groups)
            ds.(gname).groups = struct();
            for k=1:numel(ginfo.Groups)
                new_ginfo = ginfo.Groups(k);
                new_longname = ds.(gname).LongName;
                ds.(gname).groups = nc_read_file_group(filename,new_ginfo,...
                    ds.(gname).groups,new_longname);
            end
        else
            return
        end
    end
end


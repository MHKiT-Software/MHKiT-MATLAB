function ds = rotate2(data_set, out_frame)
    % Rotate a dataset to a new coordinate system.
    % 
    % Parameters
    % ----------
    % ds : xr.Dataset
    %   The dolfyn dataset (ADV or ADCP) to rotate.
    % out_frame : string {'beam', 'inst', 'earth', 'principal'}
    %   The coordinate system to rotate the data into.
    % 
    % Returns
    % -------
    % ds : xarray.Dataset
    %   The rotated dataset
    % 
    % Notes
    % -----
    % This function rotates all variables in ds.attrs.rotate_vars

    % The 'rotation chain'
    rc = {'beam', 'inst', 'earth', 'principal'};
    rot_module_dict = {'vector', 'awac', 'signature', 'ad2cp', 'rdi'};

    csin = lower(data_set.coord_sys);
    if strcmp(csin,'ship')
        csin = 'inst';
    end

    r_vec = true;
    if ~isfield(data_set,'inst2head_rotmat')
        r_vec = false;
    end
    if r_vec && ~ data_set.inst2head_rotmat_was_set
        msgtext = ['The inst2head rotation matrix exists in props, '...
            'but it was not set using `set_inst2head_rotmat.'];
        ME = MException('MATLAB:read_nortek:set_declination:rotate2'...
            ,msgtext);
        throwAsCaller(ME)
    end
    if r_vec && det(data_set.inst2head_rotmat.data) ~= 1
        msgtext = ['Invalid inst2head_rotmat (determinant != 1).'];
        ME = MException('MATLAB:read_nortek:set_declination:rotate2'...
            ,msgtext);
        throwAsCaller(ME)
    end

    if strcmp(out_frame,'principal') && ~strcmp(csin,'earth')
        warning(["You are attempting to rotate into the 'principal'"...
            " coordinate system, but the dataset is in the %s " ...
            "coordinate system. Be sure that 'principal_angle' is " ...
            "defined based on the earth coordinate system."], csin);
    end

    rmod = '';
    make_model = lower(join([data_set.attrs.inst_make,' ',...
        data_set.attrs.inst_model]));       
    for qq = 1:numel(rot_module_dict)
        ky = rot_module_dict{qq};
        if (contains(make_model,ky))
            rmod = ky;
            break;
        end
    end
    if strcmp(rmod,'')
        msgtext = ['Rotations are not defined for instrument %s',...
            make_model];
        ME = MException('MATLAB:read_nortek:set_declination:rotate2'...
            ,msgtext);
        throwAsCaller(ME)
    end
    
    % Get the 'indices' of the rotation chain
    iframe_in = find(strcmp(rc,csin));
    iframe_out = find(strcmpi(rc,out_frame));

    if isempty(iframe_in)
        msgtext = ['The coordinate system of the input dataset, %s'...
            ', is invalid.',data_set.coord_sys];
        ME = MException('MATLAB:read_nortek:set_declination:rotate2'...
            ,msgtext);
        throwAsCaller(ME)
    end
    if isempty(iframe_out)
        msgtext = ["The specified output coordinate system is "...
            "invalid, please select one of: 'beam', 'inst', "...
            "'earth', 'principal'."];
        ME = MException('MATLAB:read_nortek:set_declination:rotate2'...
            ,msgtext);
        throwAsCaller(ME)
    end

    if iframe_out == iframe_in
        fprintf('Data is already in the %s coordinate system', out_frame)
        ds = data_set;
        return;
    end

    if iframe_out > iframe_in
        reverse = false;
    else
        reverse = true;
    end

    while ~strcmpi(data_set.coord_sys, out_frame)
        inow = find(strcmp(rc,csin));
        if reverse
            func = join([rc(inow-1),'2',rc(inow)],'');
        else
            func = join([rc(inow),'2',rc(inow+1)],'');
        end
        try
            data_set = feval(func{1}, data_set, reverse, rmod);
        catch
            data_set = feval(func{1}, data_set, reverse);
        end
    end
    ds = data_set;
end

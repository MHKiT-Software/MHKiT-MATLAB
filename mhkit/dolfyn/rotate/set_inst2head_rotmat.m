function data_set = set_inst2head_rotmat(ds, rotmat)
%
%     Set the instrument to head rotation matrix for cable head Nortek ADVs.
%
%     Parameters
%     ----------
%     rotmat : float
%         3x3 rotation matrix
%
%     Returns
%     ----------
%     data_set : structure
%         Dataset with rotation matrix applied
%
%
    if ~strcmpi(ds.attrs.inst_model, 'vector')
        msgtext = ['Setting inst2head_rotmat is only supported',...
            'for Nortek Vector ADVs.'];
        ME = MException('MATLAB:dolfyn:roate:set_inst2head_rotmat'...
            ,msgtext);
        throwAsCaller(ME)
    end
    if isfield(ds,'inst2head_rotmat')
        msgtext = ['You are setting inst2head_rotmat after it has ',...
            'already been set. You can only set it once.'];
        ME = MException('MATLAB:dolfyn:roate:set_inst2head_rotmat'...
            ,msgtext);
        throwAsCaller(ME)
    end

    csin = ds.coord_sys;

    if ~strcmpi(csin, 'inst') && ~strcmpi(csin, 'beam')
        ds = rotate2(ds, 'inst');
    end

    ds.inst2head_rotmat.data = rotmat;
    ds.inst2head_rotmat.dims = {'x_star', 'x'};
    ds.inst2head_rotmat.coords.x_star = [1,2,3];
    ds.inst2head_rotmat.coords.x = [1,2,3];
    ds.attrs.inst2head_rotmat_was_set = true;

    if  ~strcmpi(csin, 'beam')
        v_shape = size(ds.vel.data);
        ds.vel.data = reshape(squeeze(ds.vel.data)*rotmat,v_shape);
    end
    if ~strcmpi(csin, 'inst') && ~strcmpi(csin, 'beam')
        ds = rotate2(ds, csin);
    end

    data_set = ds;
end


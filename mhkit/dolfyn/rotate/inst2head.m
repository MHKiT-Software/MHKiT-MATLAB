function ds = inst2head(advo, reverse)
if ~ check_inst2head_rotmat(advo)
    ds = advo;
    return
end

rotmat = advo.inst2head_rotmat.data;
if reverse
    v_shape = size(advo.vel.data);
    advo.vel.data = reshape(squeeze(advo.vel.data)*rotmat',v_shape);
else
    v_shape = size(advo.vel.data);
    advo.vel.data = reshape(squeeze(advo.vel.data)*rotmat,v_shape);
end

ds = advo;

    function result = check_inst2head_rotmat(advo)
        if ~isfield(advo, 'inst2head_rotmat')
            % This is the default value, and we do nothing.
            result = false;
            return
        end
        if ~advo.attrs.inst2head_rotmat_was_set
            msgtext = ['The inst2head rotation matrix exists in props,' ...
                'but it was not set using `set_inst2head_rotmat.'];
            ME = MException('MATLAB:dolfyn:rotate:inst2head',msgtext);
            throwAsCaller(ME)
        end
        if det(advo.inst2head_rotmat.data) ~= 1
            ME = MException('MATLAB:dolfyn:rotate:inst2head',"Invalid " + ...
                "inst2head_rotmat (determinant != 1)");
            throwAsCaller(ME)
        end
        result = true;
    end
end


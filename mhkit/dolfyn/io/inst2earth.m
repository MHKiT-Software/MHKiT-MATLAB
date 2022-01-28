function ds = inst2earth(advo,reverse)
%     Rotate data in an ADV object to the earth from the instrument
%     frame (or vice-versa).
% 
%     Parameters
%     ----------
%     advo : The adv object containing the data.
% 
%     reverse : bool (default: False)
%            If True, this function performs the inverse rotation
%            (earth->inst).
% 
%     rotate_vars : iterable
%       The list of variables to rotate. By default this is taken from
%       advo.props['rotate_vars'].
% 
%     force : Do not check which frame the data is in prior to
%       performing this rotation.

    if reverse % earth -> inst
        % The transpose of the rotation matrix gives the inverse
        % rotation, so we simply reverse the order of the einsum:
        sumstr = 'dcab,dcbe->dcae';%'jik,j...k->i...k';
        cs_now = 'earth';
        cs_new = 'inst';
    else % inst->earth
        sumstr = 'dcba,dcbe->dcae'; %'ijk,j...k->i...k'; 
        cs_now = 'inst';
        cs_new = 'earth';
    end

    if isfield(advo.attrs, 'rotate_vars')
        rotate_vars = advo.attrs.rotate_vars;
    else
        rotate_vars = {'vel'};
    end

    cs = lower(advo.coord_sys);
    if strcmp(cs,cs_new)
        return;
    elseif ~strcmp(cs, cs_now)
        msgtext = ["Data must be in the '%s' frame when using this "...
            "function.", cs_now];
        ME = MException('MATLAB:read_nortek:set_declination:rotate2'...
            ,msgtext);
        throwAsCaller(ME)
    end

    if isfield(advo, 'orientmat')
        omat = advo.orientmat.data;
    else
        if contains(lower(data_set.attrs.inst_model),'vector')
            orientation_down = advo.orientation_down;
        else
            orientation_down = nan;
        end
        omat = calc_omat(advo.heading.data,...
                         advo.pitch.data,...
                         advo.roll.data,...
                         orientation_down);
    end

    % Take the transpose of the orientation to get the inst->earth rotation
    % matrix.
    rmat = permute(omat, [1 2 4 3]);
    det_form = permute(omat, [4 3 2 1]);
    determinate = 0;
    for i = 1:length(det_form)
        determinate = determinate + det(det_form(:,:,:,i));
    end
    determinate = determinate - length(det_form);
    if determinate > 1e-3
        warning("Invalid orientation matrix (determinant != 1) " + ...
            "in inst2earth")
    end

    for qq = 1:numel(rotate_vars)
        nm = rotate_vars{qq};
        n = size(advo.(nm).data);
        n = n(3);
        if n ~= 3
            msgtext = ["The entry %s is not a vector, it cannot be "...
                "rotated.", nm];
            ME = MException('MATLAB:read_nortek:set_declination:rotate2'...
                ,msgtext);
            throwAsCaller(ME)
        end
        advo.(nm).data = tensorproduct(rmat,advo.(nm).data,sumstr);
    end
    advo.coord_sys = cs_new;
    ds = advo;
end


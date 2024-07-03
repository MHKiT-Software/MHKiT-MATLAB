function ds = earth2principal(advo,reverse)
%     Rotate data in an ADV dataset to/from principal axes. Principal
%     heading must be within the dataset.
%
%     All data in the advo.attrs['rotate_vars'] list will be
%     rotated by the principal heading, and also if the data objet has an
%     orientation matrix (orientmat) it will be rotated so that it
%     represents the orientation of the ADV in the principal
%     (reverse:earth) frame.
%
%     Parameters
%     ----------
%     advo : The adv object containing the data.
%     reverse : bool (default: False)
%            If True, this function performs the inverse rotation
%            (principal->earth).
%

% This is in degrees CW from North
ang = deg2rad(90 - advo.attrs.principal_heading);
% convert this to radians CCW from east (which is expected by
% the rest of the function)

if reverse
    cs_now = 'principal';
    cs_new = 'earth';
else
    ang = -1 * ang;
    cs_now = 'earth';
    cs_new = 'principal';
end

cs = lower(advo.coord_sys);
if strcmp(cs,cs_new)
    fprintf('Data is already in the %s coordinate system', cs)
    ds = advo;
    return;
elseif ~strcmp(cs, cs_now)
    msgtext = ["Data must be in the '%s' frame when using this "...
        "function.", cs_now];
    ME = MException('MATLAB:dolfyn:earth2principal',msgtext);
    throwAsCaller(ME)
end

ds = advo;

% Calculate the rotation matrix:
cp = cos(ang);
sp = sin(ang);
rotmat = [cp, -sp, 0; sp, cp, 0; 0, 0, 1];

% Perform the rotation:
if isfield(advo.attrs, 'rotate_vars')
    rotate_vars = advo.attrs.rotate_vars;
else
    rotate_vars = {'vel'};
end

for qq = 1:numel(rotate_vars)
    nm = rotate_vars{qq};
    shape = size(advo.(nm).data);
    l = length(shape);
    otherdims = repmat({':'},1,l-1);
    if l == 4
        sumstr = 'da,cba->cbd';
    elseif l == 3
        sumstr = 'da,ca->cd';
    end
    rot = tensorproduct(rotmat(1:2,1:2),...
        squeeze(advo.(nm).data(otherdims{:},1:2)),sumstr);
    shape(end) = 2;
    ds.(nm).data(otherdims{:},1:2) = reshape(rot,shape);
end

ds = set_coords(ds, cs_new);
ds.coord_sys = cs_new;
ds.attrs.coord_sys = cs_new;

end


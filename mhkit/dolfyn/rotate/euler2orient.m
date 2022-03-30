function omat = euler2orient(heading, pitch, roll)
    % For Nortek data only.
    % The heading, pitch, roll used here are from the Nortek binary files.
    %
    % Heading input is clockwise from North
    % Returns a rotation matrix that rotates earth (ENU) -> inst.
    % This is based on the Nortek `Transforms.m` file, available in
    % the refs folder.
    heading = deg2rad(heading);
    roll = deg2rad(roll);
    pitch = deg2rad(pitch);

    % The definition of heading below is consistent with the right-hand-rule;
    % heading is the angle positive counterclockwise from North of the y-axis.
    %
    % This also involved swapping the sign on sh in the def of omat
    % below from the values provided in the Nortek Matlab script.
    heading = (pi/2 - heading);

    ch = cos(heading);
    sh = sin(heading);
    cp = cos(pitch);
    sp = sin(pitch);
    cr = cos(roll);
    sr = sin(roll);

    % Note that I've transposed these values (from what is defined in
    % Nortek matlab script), so that the omat is earth->inst
    omat = zeros([length(sh),1,3,3], 'double');
    omat(:,:,1,1) = ch .* cp;
    omat(:,:,1,2) = -ch .* sp .* sr - sh .* cr;
    omat(:,:,1,3) = -ch .* cr .* sp + sh .* sr;
    omat(:,:,2,1) = sh .* cp;
    omat(:,:,2,2) = -sh .* sp .* sr + ch .* cr;
    omat(:,:,2,3) = -sh .* cr .* sp - ch .* sr;
    omat(:,:,3,1) = sp;
    omat(:,:,3,2) = sr .* cp;
    omat(:,:,3,3) = cp .* cr;
end


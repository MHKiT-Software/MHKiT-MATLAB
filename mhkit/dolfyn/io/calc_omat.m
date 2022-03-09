function omat = calc_omat(hh, pp, rr, orientation_down)
    hh_check = isnan(hh);
    pp_check = isnan(pp);
    rr_check = isnan(rr);

    if hh_check(end) && pp_check(end) && rr_check(end)
        % The end of the data may not have valid orientations
        last_non_NaN_index_hh = find(~isnan(hh), 1, 'last');
        last_non_NaN_index_pp = find(~isnan(pp), 1, 'last');
        last_non_NaN_index_rr = find(~isnan(rr), 1, 'last');
        lastgd = min([last_non_NaN_index_hh,...
                     last_non_NaN_index_pp,...
                     last_non_NaN_index_rr]);
        hh(lastgd:end) = hh(lastgd);
        pp(lastgd:end) = pp(lastgd);
        rr(lastgd:end) = rr(lastgd);
    end
    if ~isnan(orientation_down)
        % For Nortek Vector ADVs: 'down' configuration means the head was
        % pointing 'up', where the 'up' orientation corresponds to the
        % communication cable being up.  
        rr(orientation_down) = rr(orientation_down) + 180;
    end

    % For Nortek data only.
    % The heading, pitch, roll used here are from the Nortek binary files.
    %
    % Heading input is clockwise from North
    % Returns a rotation matrix that rotates earth (ENU) -> inst.
    % This is based on the Nortek `Transforms.m` file, available in
    % the refs folder.
    heading = deg2rad(hh);
    roll = deg2rad(rr);
    pitch = deg2rad(pp);

    % The definition of heading below is consistent with the 
    % right-hand-rule; heading is the angle positive counterclockwise
    % from North of the y-axis.
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

% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


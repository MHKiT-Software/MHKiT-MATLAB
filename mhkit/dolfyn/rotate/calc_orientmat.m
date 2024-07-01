function rotmat = calc_orientmat(adcpo)
%%%%%%%%%%%%%%%%%%%%
<<<<<<< HEAD
%    Calculate the orientation matrix using the raw 
%    heading, pitch, roll values from the RDI binary file.
%     
=======
%    Calculate the orientation matrix using the raw
%    heading, pitch, roll values from the RDI binary file.
%
>>>>>>> master
% Parameters
% ------------
%    adcpo : The ADP object containing the data.
%
%    ## RDI-ADCP-MANUAL (Jan 08, section 5.6 page 18)
%    The internal tilt sensors do not measure exactly the same
%    pitch as a set of gimbals would (the roll is the same). Only in
%    the case of the internal pitch sensor being selected (EZxxx1xxx),
%    the measured pitch is modified using the following algorithm.
%
%        P = arctan[tan(Tilt1)*cos(Tilt2)]    (Equation 18)
%
%    Where: Tilt1 is the measured pitch from the internal sensor, and
%    Tilt2 is the measured roll from the internal sensor The raw pitch
%    (Tilt 1) is recorded in the variable leader. P is set to 0 if the
%    "use tilt" bit of the EX command is not set."""
%
% Returns
% ---------
<<<<<<< HEAD
%     rotmat: array 
%         orientation matrix
%        
=======
%     rotmat: array
%         orientation matrix
%
>>>>>>> master
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    r = deg2rad(adcpo.roll.data);
    p = atan(tan(deg2rad(adcpo.pitch.data)).* cos(r));
    h = deg2rad(adcpo.heading.data);

    if contains(lower(adcpo.attrs.inst_make), 'rdi')
        if strcmp(adcpo.attrs.orientation,'up')
            % RDI-ADCP-MANUAL (Jan 08, section 5.6 page 18)
            % Since the roll describes the ship axes rather than the
            % instrument axes, in the case of upward-looking
            % orientation, 180 degrees must be added to the measured
            % roll before it is used to calculate M. This is equivalent
            % to negating the first and third columns of M. R is set
            % to 0 if the "use tilt" bit of the EX command is not set.
            r = r + pi;
        end

        if strcmp(adcpo.attrs.coord_sys, 'ship') &&...
                adcpo.attrs.use_pitchroll
            r(:) = 0;
            p(:) = 0;
        end
    end

    ch = cos(h);
    sh = sin(h);
    cr = cos(r);
    sr = sin(r);
    cp = cos(p);
    sp = sin(p);

    rotmat = zeros(length(r),1,3,3);
    rotmat(:,:,1,1) = ch .* cr + sh .* sp .* sr;
    rotmat(:,:,2,1) = sh .* cp;
    rotmat(:,:,3,1) = ch .* sr - sh .* sp .* cr;
    rotmat(:,:,1,2) = -sh .* cr + ch .* sp .* sr;
    rotmat(:,:,2,2) = ch .* cp;
    rotmat(:,:,3,2) = -sh .* sr - ch .* sp .* cr;
    rotmat(:,:,1,3) = -cp .* sr;
    rotmat(:,:,2,3) = sp;
<<<<<<< HEAD
    rotmat(:,:,3,3) = cp .* cr;   
=======
    rotmat(:,:,3,3) = cp .* cr;
>>>>>>> master

    rotmat = permute(rotmat, [1 2 4 3]);
end


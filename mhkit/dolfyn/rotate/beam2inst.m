function ds = beam2inst(dat,reverse,rmod)
% Rotate velocities from beam to instrument coordinates.
%
%    Parameters
%    ----------
%    dat : The adv object containing the data.
%        The ADCP dataset
%    reverse : bool 
%        If True, this function performs the inverse rotation (inst->beam).
%    rmod: string
%        awac, vector, signature, rdi
%__________________________________________________________________________

% There are two paths depending on the make of the data
if strcmpi(rmod,'vector')
    % Order of rotations matters
    % beam->head(ADV instrument head)->inst(ADV battery case|imu)
    if reverse
        % First rotate velocities from ADV inst frame back to head frame
        dat = inst2head(dat, reverse);
        % Now rotate from the head frame to the beam frame
        ds = rotate(dat, reverse);
    
    % inst(ADV battery case|imu)->head(ADV instrument head)->beam
    else
        % First rotate velocities from beam to ADV head frame
        dat = rotate(dat, reverse);
        % Then rotate from ADV head frame to ADV inst frame
        ds = inst2head(dat, reverse);
    end
else
    ds = rotate(dat, reverse);
end
%__________________________________________________________________________

    function data = rotate(advo,reverse)

        if reverse && ~strcmpi(advo.coord_sys, 'inst')
            msgtext = ['The input must be in inst coordinates.'];
            ME = MException('MATLAB:dolfyn:rotate:beam2inst',msgtext);
            throwAsCaller(ME)
        end
        if ~reverse && ~strcmpi(advo.coord_sys, 'beam')
            msgtext = ['The input must be in beam coordinates.'];
            ME = MException('MATLAB:dolfyn:rotate:beam2inst',msgtext);
            throwAsCaller(ME)
        end
        
        if ~isfield(advo, 'beam2inst_orientmat')
            msgtext = ['Data set is missing the beam2inst orientation ' ...
                'matrix'];
            ME = MException('MATLAB:dolfyn:rotate:beam2inst',msgtext);
            throwAsCaller(ME)
        end
        rotmat = advo.beam2inst_orientmat.data;
        
        cs = 'inst';
        if reverse
            % Can't use transpose because rotation is not between
            % orthogonal coordinate systems
            rotmat = inv(rotmat);
            cs = 'beam';
        end
        
        for qq = 1:numel(advo.attrs.rotate_vars)
            nm = advo.attrs.rotate_vars{qq};
            if startsWith(nm, 'vel')
                if length(size(advo.(nm).data)) == 3
                    advo.(nm).data = tensorproduct(...
                        rotmat,advo.(nm).data, 'ab,cda->cdb');
                else
                    advo.(nm).data = tensorproduct(...
                        rotmat,advo.(nm).data, 'ab,cdea->cdeb');
                end
            end
        end
        
        data = set_coords(advo, cs);
    end

end


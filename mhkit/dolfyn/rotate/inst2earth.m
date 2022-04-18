function ds = inst2earth(advo,reverse,make)
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
%     make : char
%       vector, signature, rdi

if strcmpi(make,'awac') 
    ds = inst2earth_awac();
elseif strcmpi(make,'vector')
    ds = inst2earth_vector();
elseif strcmpi(make, 'signature')
    ds = inst2earth_signature();
elseif strcmpi(make, 'rdi')
    ds = inst2earth_rdi();
end

    function out = inst2earth_awac()
        if reverse % earth -> inst
            % The transpose of the rotation matrix gives the inverse
            % rotation, so we simply reverse the order of the einsum:
            sumstr = 'dcab,dceb->dcea';%'jik,j...k->i...k';
            cs_now = 'earth';
            cs_new = 'inst';
        else % inst->earth
            sumstr = 'dcba,dceb->dcea'; %'ijk,j...k->i...k'; 
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
            ME = MException('MATLAB:dolfyn:inst2earth',msgtext);
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
            n = n(end);
            if n ~= 3
                msgtext = ["The entry %s is not a vector, it cannot be "...
                    "rotated.", nm];
                ME = MException('MATLAB:inst2earth',msgtext);
                throwAsCaller(ME)
            end
            advo.(nm).data = tensorproduct(rmat,advo.(nm).data,sumstr);
        end
        advo.coord_sys = cs_new;
        advo = set_coords(advo, cs_new);
        out = advo;
    end

    function out = inst2earth_vector()
        if reverse % earth -> inst
            % The transpose of the rotation matrix gives the inverse
            % rotation, so we simply reverse the order of the einsum:
            sumstr = 'dab,dcb->dca';%'jik,j...k->i...k';
            cs_now = 'earth';
            cs_new = 'inst';
        else % inst->earth
            sumstr = 'dba,dcb->dca'; %'ijk,j...k->i...k'; 
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
            ME = MException('MATLAB:dolfyn:inst2earth',msgtext);
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
        rmat = squeeze(rmat);
    
        for qq = 1:numel(rotate_vars)
            nm = rotate_vars{qq};
            n = size(advo.(nm).data);
            n = n(end);
            if n ~= 3
                msgtext = ["The entry %s is not a vector, it cannot be "...
                    "rotated.", nm];
                ME = MException('MATLAB:inst2earth',msgtext);
                throwAsCaller(ME)
            end
            advo.(nm).data = tensorproduct(rmat,advo.(nm).data,sumstr);
        end
        advo.coord_sys = cs_new;
        advo = set_coords(advo, cs_new);
        out = advo;
    end

    function out = inst2earth_signature()
        if reverse % earth -> inst
            % The transpose of the rotation matrix gives the inverse
            % rotation, so we simply reverse the order of the einsum:
            sumstr = 'dac,dbc->dba';%'jik,j...k->i...k';
            cs_now = 'earth';
            cs_new = 'inst';
        else % inst->earth
            sumstr = 'dca,dbc->dba'; %'ijk,j...k->i...k'; 
            cs_now = 'inst';
            cs_new = 'earth';
        end

        if isfield(advo.attrs, 'orientation')
            if strcmpi(advo.attrs.orientation, 'down')
                down = true;
            else
                down = false;
            end
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
            omat = euler2orient( advo.heading.data,...
                                 advo.pitch.data,...
                                 advo.roll.data);
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
    
        % The dictionary of rotation matrices for different sized arrays.
        rmd_3 = rmat;
        
        % The 4-row rotation matrix assume that rows 0,1 are u,v,
        % and 2,3 are independent estimates of w.
        tmp = zeros([length(rmat),1,4,4],"double");
        tmp(:,:,1:3,1:3) = rmd_3;
        % Copy row 3 to row 4
        tmp(:,:,1:2,4) = rmat(:,:,1:2,3);
        tmp(:,:,4,4) = rmat(:,:,3,3);
        % Extend rows 1 and 2
        tmp(:,:,3,1) = rmat(:,:,3,1)./2.;
        tmp(:,:,4,1) = rmat(:,:,3,1)./2.;

        tmp(:,:,3,2) = rmat(:,:,3,2)./2.;
        tmp(:,:,4,2) = rmat(:,:,3,2)./2.;

        rmd_3 = squeeze(rmd_3);
        rmd_4 = tmp;

        if reverse
            % 3-element inverse handled by sumstr definition (transpose)
            step1 = squeeze(permute(rmd_4,[3,4,1,2]));
            % matlab cant invert n-d matricies 
            shape = size(step1);
            step2 = zeros(shape);
            for kk = 1:shape(3)
                step2(:,:,kk) = inv(step1(:,:,kk));
            end
            step3 = reshape(step2,[shape(1:2),1,shape(3)]);
            rmd_4 = permute(step3,[4,3,1,2]);
        end

        for qq = 1:numel(advo.attrs.rotate_vars)
            nm = advo.attrs.rotate_vars{qq};
            dat = advo.(nm).data;
            shape = size(dat);
            n = shape(end);
            % Nortek documents sign change for upside-down instruments
            if down                
                if ~reverse
                    for kk = 2:n
                        dat(:,:,:,kk) = dat(:,:,:,kk) .* -1.;
                    end
                    if n == 3
                        dat = tensorproduct(rmd_3,dat,sumstr);
                    elseif n == 4
                        dat = tensorproduct(rmd_4,dat,'dcba,dceb->dcea');
                    else
                        msgtext = ["The entry %s is not a vector, it " +...
                            "cannot be rotated.", nm];
                        ME = MException('MATLAB:inst2earth',msgtext);
                        throwAsCaller(ME)
                    end
                else
                    if n == 3
                        dat = tensorproduct(rmd_3,dat,sumstr); 
                    elseif n == 4
                        dat = tensorproduct(rmd_4,dat,'dcba,dceb->dcea');                        
                    else
                        msgtext = ["The entry %s is not a vector, it " +...
                            "cannot be rotated.", nm];
                        ME = MException('MATLAB:inst2earth',msgtext);
                        throwAsCaller(ME)
                    end
                    for kk = 2:n
                        dat(:,:,:,kk) = dat(:,:,:,kk) .* -1.;
                    end
                end
            else % 'up' and AHRS
                if n == 3                    
                    dat = tensorproduct(rmd_3,dat,sumstr);
                elseif n == 4
                    dat = tensorproduct(rmd_4,dat,'dcba,dceb->dcea');
                else
                    msgtext = ["The entry %s is not a vector, it " + ...
                            "cannot be rotated.", nm];
                        ME = MException('MATLAB:inst2earth',msgtext);
                        throwAsCaller(ME)
                end
            end 
            advo.(nm).data = dat;
        end

        advo = set_coords(advo, cs_new);
        advo.coord_sys = cs_new;
        out = advo;
    end

    function out = inst2earth_rdi()
        % Rotate velocities from the instrument to earth coordinates.
        % 
        % This function also rotates data from the 'ship' frame, into the
        % earth frame when it is in the ship frame (and
        % ``adcpo.use_pitchroll == 'yes'``). It does not support the
        % 'reverse' rotation back into the ship frame.
        % 
        % Parameters
        % ----------
        % adpo : The ADP object containing the data.
        % 
        % reverse : bool (default: False)
        %        If True, this function performs the inverse rotation
        %        (earth->inst).
        % fixed_orientation : bool (default: False)
        %     When true, take the average orientation and apply it over the
        %     whole record.
        % force : bool (default: False)
        %     When true do not check which coordinate system the data is in
        %     prior to performing this rotation.
        % 
        % Notes
        % -----
        % The rotation matrix is taken from the Teledyne RDI ADCP Coordinate
        % Transformation manual January 2008
        csin = lower(advo.coord_sys);
        cs_allowed = {'inst', 'ship'};
        if reverse
            cs_allowed = {'earth'};
        end

        if ~any(strcmpi(csin, cs_allowed))
            msgtext = ["Invalid rotation for data in %s coordinate system."...
                , csin];
            ME = MException('MATLAB:read_nortek:set_declination:rotate2'...
                ,msgtext);
            throwAsCaller(ME)
        end

        if isfield(advo, 'orientmat')
            omat = advo.orientmat.data;
        else
            omat = calc_orientmat(advo);
        end

        % rollaxis gives transpose of orientation matrix.
        % The 'rotation matrix' is the transpose of the 'orientation matrix'
        % NOTE: the double 'rollaxis' within this function, and here, has
        % minimal computational impact because np.rollaxis returns a
        % view (not a new array)
        rmat = permute(omat, [1 2 4 3]);
        if reverse % earth -> inst
            sumstr = 'dcab,dceb->dcea';%'jik,j...k->i...k';
            cs_new = 'inst';
        else % inst->earth
            sumstr = 'dcba,dceb->dcea'; %'ijk,j...k->i...k';
            cs_new = 'earth';
        end

        % Only operate on the first 3-components, b/c the 4th is err_vel
        for qq = 1:numel(advo.attrs.rotate_vars)
            nm = advo.attrs.rotate_vars{qq};
            dat = advo.(nm).data(:,:,:,1:3);
            advo.(nm).data(:,:,:,1:3) = tensorproduct(rmat,dat,sumstr);            
        end

        advo = set_coords(advo, cs_new);
        advo.coord_sys = cs_new;
        out = advo;
    end
end


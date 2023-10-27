function out = set_coords(ds, ref_frame)
%     Checks the current reference frame and adjusts xarray coords/dims 
%     as necessary. Makes sure assigned dataarray coordinates match what
%     DOLfYN is reading in.
    make = join([ds.attrs.inst_make," ", ds.attrs.inst_model]);
    XYZ = {'X', 'Y', 'Z'};
    ENU = {'E', 'N', 'U'};
    beam = 1:min(size(squeeze(ds.vel.data)));
    principal = {'streamwise', 'x-stream', 'vert'};

    if contains(lower(make), 'rdi')
        inst  = {'X', 'Y', 'Z', 'err'};
        earth = {'E', 'N', 'U', 'err'};
        princ = {'streamwise', 'x-stream', 'vert', 'err'};

    elseif contains(lower(make), 'nortek')
        if contains(lower(make), 'signature') || ...
                contains(lower(make), 'ad2cp')
            inst  = {'X', 'Y', 'Z1', 'Z2'};
            earth = {'E', 'N', 'U1', 'U2'};
            princ = {'streamwise', 'x-stream', 'vert1', 'vert2'};
        else
            % AWAC or Vector
            inst = XYZ;
            earth = ENU;
            princ = principal;
        end
    end

    %  update 'orient' and 'orientIMU' dimensions
    if strcmpi(ref_frame, 'beam')
        ds.coords.dir = beam;
    elseif strcmpi(ref_frame, 'inst')
        ds.coords.dir = inst;
    elseif strcmpi(ref_frame, 'ship')
        ds.coords.dir = inst;
    elseif strcmpi(ref_frame, 'earth')
        ds.coords.dir = earth;
    elseif strcmpi(ref_frame, 'principal')
        ds.coords.dir = princ;
    end

    if isfield(ds.coords, 'dirIMU')
        if strcmpi(ref_frame, 'earth')
            ds.coords.dirIMU = ENU;
        elseif strcmpi(ref_frame, 'principal')
            ds.coords.dirIMU = principal;
        else
            ds.coords.dirIMU = XYZ;
        end
    end
    ds.attrs.coord_sys = ref_frame;
    ds.coord_sys = ref_frame;

    % In xarray this is adequate. When you change ds.coords.dir it
    % propogates that change to each field that has dir as a coordinate.
    % Here, however, we must manaully propogate that change
    fn = fieldnames(ds); 
    for qq = 1:numel(fn)
        field = fn{qq};
        if isfield(ds.(field),'coords')
            % This is a data field
            if isfield(ds.(field).coords, 'dir')
                % replace dir
                ds.(field).coords.dir = ds.coords.dir;
            end
            if isfield(ds.(field).coords, 'dirIMU')
                % replace dirIMU
                ds.(field).coords.dirIMU = ds.coords.dirIMU;
            end
        end
    end
         
    tag = {'', '_echo', '_bt'};
    for qq = 1:3
        if isfield(ds.attrs, join(["coord_sys_axes",tag{qq}],''))
            ds.attrs = rmfield(ds.attrs,...
                join(["coord_sys_axes",tag{qq}],''));
        end
    end

    out = ds;
end


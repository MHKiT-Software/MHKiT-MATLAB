function out = nan_beyond_surface(ds, options)
% Mask the values of 3D data (vel, amp, corr, echo) that are beyond
% the surface.
%
% Parameters
% ----------
% ds : Dataset
%   The adcp dataset to clean
% val : nan or numeric
%   Specifies the value to set the bad values to (default np.nan).
%
% Returns
% -------
% ds : Dataset
%   The adcp dataset where relevant arrays with values greater than
%   `depth` are set to NaN
%
% Notes
% -----
% Surface interference expected to happen at
% distance > range * cos(beam angle) - cell size

    arguments
        ds
        options.val = nan;
    end

    fn = fieldnames(ds);
    r = {};
    for qq = 1:numel(fn)
        field = fn{qq};
        if isfield(ds.(field),'coords')
            % This is a data field
            for kk = 1:numel(ds.(field).dims)
                % Check the coords for any range values
                if contains(ds.(field).dims{kk},'range')
                    r{end+1} = field;
                end
            end
        end
    end

    if contains(lower(ds.attrs.inst_make),'nortek')
        beam_angle = 25 * (pi/180);
    else  % TRDI
        try
            beam_angle = ds.attrs.beam_angle;
        catch
            beam_angle = 20 * (pi/180);
        end
    end

    % Surface interference distance calculated from distance of
    % transducers to surface
    if isfield(ds.attrs, 'h_deploy')
        range_limit = ((ds.depth.data-ds.attrs.h_deploy).*cos(beam_angle) ...
                        - ds.attrs.cell_size) + ds.attrs.h_deploy;
    else
        range_limit = ds.depth.data .* cos(beam_angle) - ...
            ds.attrs.cell_size;
    end

    % Create a 2d matrix of to compare each range against the range limit
    ds_coords_range_rep = repmat(ds.coords.range', numel(range_limit), 1);

    % Create a 2d matrix of binary values
    % 1: Below the surface interference distance
    % 0: Above the surface interference distance
    bds = (ds_coords_range_rep > range_limit);


    % Echosounder data gets trimmed at water surface
    if any(strcmp(r,'echo'))
        bds_echo = ds.coords.range_echo > ds.depth;
        ech_shape = size(bds_echo);
        bds_echo = reshape(bds_echo,[ech_shape(1),1,ech_shape(2)]);
        ds.echo.data(ech_shape) = options.val;
        r(strncmpi(r,'echo',1)) = [];
    end

    % Correct rest of "range" data for surface interference
    for qq = 1:numel(r)
        a = ds.(r{qq}).data;
        if length(size(a)) > 3
            for kk = 1:size(a,4)
                sub_list = a(:,:,:,kk);
                try  % float dtype
                    sub_list(bds) = options.val;
                catch  % int dtype
                    sub_list(bds) = 0;
                end
                a(:,:,:,kk) = sub_list;
            end
        else
            try  % float dtype
                a(bds) = options.val;
            catch  % int dtype
                a(bds) = 0;
            end
        end
        ds.(r{qq}).data = a;
    end

    out = ds;

end

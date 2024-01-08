function out = correlation_filter(ds, options)
% Filters out velocity data where correlation is below a
% threshold in the beam correlation data.
%
% Parameters
% ----------
% ds : Dataset
%   The adcp dataset to clean.
% thresh : numeric
%   The maximum value of correlation to screen, in counts or %, default is
%   50
% val : numeric
%   Value to set masked correlation data to, default is nan
%
% Returns
% -------
% ds : xarray.Dataset
%  Velocity data with low correlation values set to `val`
%
% Notes
% -----
% Does not edit correlation or amplitude data.

    arguments
        ds;
        options.thresh = 50;
        options.val = nan;
    end

    % 4 or 5 beam
    if isfield(ds,'vel_b5')
        tag = {'', '_b5'};
    else
        tag = {''};
    end

    % copy original ref frame
    coord_sys_orig = ds.coord_sys;

    % correlation is always in beam coordinates
    ds = rotate2(ds, 'beam');
    for qq = 1:numel(tag)
        tg = string(tag{qq});
        field = strcat("corr",tg);
        mask = ds.(field).data < options.thresh;
        field = strcat("vel",tg);
        ds.(field).data(mask) = options.val;
        comment = "Filtered of data with a correlation value below " + ...
            string(options.thresh) + string(ds.corr.units);
        ds.(field).comment = comment;
    end

    out = rotate2(ds, coord_sys_orig);

end


function out = find_surface(ds, options)
%TODO: update comments here since they are just the python comments
%unchanged
% Find the surface (water level or seafloor) from amplitude data and
% adds the variable "depth" to the input Dataset.
%
% Parameters
% ----------
% vds : struct (from a fieldname in dataset in create_dataset function)
%   The full adcp dataset
% options.thresh : int
%   Specifies the threshold used in detecting the surface. Default = 10
%     (The amount that amplitude must increase by near the surface for it to
%     be considered a surface hit)
% options.nfilt : int
%   Specifies the width of the median filter applied, must be odd.
%     Default is None
%
% Returns
% -------
% None, operates "in place"

    % set default arguments
    arguments
        ds;
        options.thresh {mustBeInteger} = 10;
        options.nfilt {mustBeInteger} = 0;
        %TODO check nfilt is odd
    end

    %TODO/NOTES: So far, this is the best attempt at translating this
    %function from python into matlab...Since python arrays are 0 indexed
    %and matlab arrays are 1 indexed some of the indexing in the
    %translation gets a bit messy. A lot of those indexes need to be double
    %checked, and there are likely a lot of off-by-one type errors.
    %Further, this has en error in it at "d1 = ds.coords.range(inds)" that
    %makes me think some of the indecies need to change even more than 1?
    %In matlab, the data has dimension Lx1x28xM and I think the second
    %dimension, 1, gets squeezed out of the python data. In other words,
    %the line below should read "[~, inds] = max(ds.amp.data, [], 3,
    %"linear")" to skip over the 1 index. (Not sure on this, but just my
    %initial thoughts on the first error in this function)

    % This finds the indices of the maximum of the echo profile:
    [~, inds] = max(ds.amp.data, [], 2, "linear")
    % This finds the first point that increases (away from the profiler) in
    % the echo profile
    % almost the same as python call, just use axis 2 since matlab 1
    % indexed.
    edf = diff(cast(ds.amp.data,"int16"), 2);
    endint = size(ds.vel.data,3);
    % use matlab : operator to get close to np.arange
    inds2 = max(uint8(edf<0) .* reshape(uint8(1:endint), 1, [], 1), [], 2);
    % Calculate the depth of these quantities
    d1 = ds.coords.range(inds);
    d2 = ds.coords.range(inds2);
    % Combine them:
    D = vertcat(d1, d2);
    % Take the median value as the estimate of the surface:
    d = median(D, 1);

    % Throw out values that do not increase near the surface by *thresh*
    for ip = 1:ds.vel.dims(2)
        itmp = min(inds(:, ip));
        if all(edf(itmp:end, :, ip) < thresh)
            d(ip) = NaN;
        end
    end

    if exists('nfilt', 'var') && nfilt ~= 0
        dfilt = medfilt2(d, [1 nfilt]);
        dfilt(conv2(isnan(d), ones(1, nfilt) / nfilt, 'same') > 4) = NaN;
        dfilt(dfilt == 0) = NaN;
        d = dfilt;
    end

    ds.depth.data = d;
    ds.depth.dims = {'time'};
    ds.depth.units='m';

    out = ds;
end


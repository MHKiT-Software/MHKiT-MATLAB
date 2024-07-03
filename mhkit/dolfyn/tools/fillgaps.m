function out = fillgaps(a, extrapFlg)
% Linearly fill NaN value in an array.
%
% Parameters
% ----------
% a : array
%   The array to be filled.
% extrapFlg : bool
%   Whether to extrapolate if NaNs are found at the ends of the
%   array.
%
% Notes
% -----
% This function interpolates assuming spacing/timestep between
% successive points is constant. If the spacing is not constant, use
% interpgaps.
%
    gd = find(~isnan(a));

    % Extrapolate if requested
    if extrapFlg && ~isempty(gd)
        if gd(1) ~= 1
            a(1:gd(1)) = a(gd(1));
        end
        if gd(end) ~= length(a)
            a(gd(end):end) = a(gd(end));
        end
    end
    % Main loop
    if length(gd) > 1
        inds = find(bitand((1 < diff(gd)),(diff(gd)<=inf)));
        for i2 = 1:length(inds)
            ii = gd(inds(i2))+1:gd(inds(i2)+1)-1;
            a(ii) = diff([a(gd(inds(i2))),a(gd(inds(i2)+1))]) *...
                (1:length(ii)) / (length(ii)+1) + a(gd(inds(i2)));
        end
    end
    out = a;
end

% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


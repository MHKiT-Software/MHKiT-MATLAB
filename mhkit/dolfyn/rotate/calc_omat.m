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

    omat = euler2orient(hh, pp, rr);
end

% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


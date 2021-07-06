function region = region_selection(lat_lon)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     Returns the name of the predefined region in which the given coordinates reside.
%     Can be used to check if the passed lat/lon pair is within the WPTO hindcast dataset. 
%     
%     Parameters
%     ----------
%         lat_lon : matrix
%             Latitude (first column) and longitude (second column) coordinates as numerics
%     
%     Returns
%     -------
%         region : string
%             Name of predefined region for given coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~isa(lat_lon,'numeric')
        error("ERROR: lat_lon must be a numeric matrix")
    end

    regions.Hawaii.lat = [15.0,15.0,27.000002,27.000002]';
    regions.Hawaii.lon = [-164.0,-151.0,-151.0,-164.0]';
    regions.West_Coast.lat = [30.0906,30.0906, 48.8641,48.8641]';
    regions.West_Coast.lon = [-130.072,-116.899,-116.899,-130.072]';
    regions.Atlantic.lat = [24.382,24.382, 44.8247,44.8247]';
    regions.Atlantic.lon = [-81.552,-65.721, -65.721,-81.552]';

    fns = fieldnames(regions);
    for i = 1:length(fns)
        in = inpolygon(lat_lon(:,1),lat_lon(:,2),regions.(fns{i}).lat,regions.(fns{i}).lon);
        if mean(in) == 1
            region = fns{i};
        end
    end
    if exist('region','var')==0
        error("WARNING: coordinates out of bounds OR not within the same region!")
    end
end
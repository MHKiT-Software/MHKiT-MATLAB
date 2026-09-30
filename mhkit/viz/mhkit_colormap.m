function map = mhkit_colormap(name, n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Return an N x 3 RGB colormap from a name, like the built-in colormap
% function, but covering MATLAB, cmocean, and viridis colormaps
%
% Names are case insensitive. Flipped colormaps use the conventions of
% their sources: cmocean names accept the MATLAB cmocean
% "-<colormap name>" form, e.g. "-thermal", and cmocean names and
% "viridis" accept the Python "<colormap name>_r" form, e.g.
% "thermal_r" or "viridis_r". To flip a MATLAB colormap pass the
% matrix, e.g. flipud(parula(256)).
%
% cmocean colormaps: https://matplotlib.org/cmocean/
% MATLAB colormaps: https://www.mathworks.com/help/matlab/ref/colormap.html
%
% Parameters
% ------------
%     name : string or matrix (optional)
%         Colormap name, one of:
%           "viridis"
%           any cmocean name: "thermal", "haline", "solar", "ice", "gray",
%             "oxy", "deep", "dense", "algae", "matter", "turbid", "speed",
%             "amp", "tempo", "rain", "phase", "topo", "balance", "delta",
%             "curl", "diff", "tarn"
%           any MATLAB colormap name: "parula", "turbo", "jet", "hot",
%             "cool", "bone", "gray", "sky", "abyss", etc.
%         An N x 3 RGB matrix is returned unchanged, so either form can be
%         passed through. Default "viridis"
%     n : integer (optional)
%         Number of colors. Default 256
%
% Returns
% ---------
%     map : matrix
%         n x 3 RGB colormap with values in [0, 1]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    name = "viridis"
    n (1,1) {mustBeInteger, mustBePositive} = 256
end

if isnumeric(name)
    if size(name, 2) ~= 3 || any(name(:) < 0) || any(name(:) > 1)
        error('MHKiT:mhkit_colormap:InvalidInput', ...
            'A colormap matrix must be N x 3 with values in [0, 1]');
    end
    map = name;
    return
end

if ~(ischar(name) || (isstring(name) && isscalar(name)))
    error('MHKiT:mhkit_colormap:InvalidInput', ...
        'name must be a colormap name or an N x 3 RGB matrix');
end
name = lower(strtrim(char(name)));
% MATLAB cmocean flips with a leading "-", Python cmocean and matplotlib
% with a trailing "_r"
minus_form = startsWith(name, '-');
r_form = endsWith(name, '_r');
flip = minus_form || r_form;
base_name = regexprep(name, '^-|_r$', '');

cmocean_names = {'thermal', 'haline', 'solar', 'ice', 'gray', 'oxy', 'deep', ...
    'dense', 'algae', 'matter', 'turbid', 'speed', 'amp', 'tempo', 'rain', ...
    'phase', 'topo', 'balance', 'delta', 'curl', 'diff', 'tarn'};

if ismember(base_name, cmocean_names)
    if flip
        map = cmocean(['-' base_name], n);
    else
        map = cmocean(base_name, n);
    end
elseif strcmp(base_name, 'viridis') && ~minus_form
    map = viridis_colormap(n);
    if flip
        map = flipud(map);
    end
elseif flip
    error('MHKiT:mhkit_colormap:InvalidInput', ...
        ['"-<colormap name>" flips cmocean colormaps and "<colormap name>_r" flips ' ...
         'cmocean and viridis colormaps. Flip "%s" by passing the matrix, e.g. flipud(%s(256))'], ...
        base_name, base_name);
elseif exist(name, 'file') == 2 || exist(name, 'builtin') == 5
    % MATLAB colormap functions all take the number of colors
    map = feval(name, n);
    if ~isnumeric(map) || size(map, 2) ~= 3
        error('MHKiT:mhkit_colormap:InvalidInput', '"%s" is not a colormap function', name);
    end
else
    error('MHKiT:mhkit_colormap:InvalidInput', ...
        'Unknown colormap "%s". Use "viridis", a cmocean name, or a MATLAB colormap name', name);
end

end

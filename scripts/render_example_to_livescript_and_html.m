function render_example_to_livescript_and_html(names, opts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Rebuild example .mlx and .html files from their .m source
%
% The .m file is the source of truth: edit it, then run this function to
% regenerate the other two formats. For each example this converts
% <name>.m to <name>.mlx, preserving publish markup such as %% headings,
% |monospace|, *bold*, $latex$, and <<image.png>>, executes the live
% script so its outputs and figures are stored in the .mlx, and exports
% <name>.html, running the script a second time so the HTML figures are
% full width. Scripts run from the examples folder, so relative data paths
% behave as they do in the Live Editor. Examples that need Python,
% Simulink, or network access must be built on a machine that has them.
%
% Figures use the light theme regardless of the OS or MATLAB theme, so
% committed output looks the same for everyone, and are rendered at
% FigureScale times the normal pixel density so they stay sharp on high
% DPI displays. The scale must be set before any Live Editor work happens
% in the MATLAB session, so run this from a fresh session, for example
% from the repository root:
%
%     matlab -batch "addpath('scripts'); render_example_to_livescript_and_html('strain_measurement_example')"
%
% Without a desktop (matlab -batch) the figure snapshots stored in the
% .mlx are limited to 469 px wide (times FigureScale). The HTML is not
% affected. For full width figures in the .mlx, open it in the MATLAB
% desktop, run it, and save.
%
% The conversion and execution steps use undocumented functions in the
% matlab.internal.liveeditor package. They are what the Live Editor uses
% internally and have been stable for many releases, but MathWorks does
% not guarantee them.
%
% Parameters
% ------------
% names : string, char, or cell array of char (optional)
%   Example names, with or without the .m extension, for example
%   'strain_measurement_example' or {'adv_example', 'adcp_example'}.
%   Default {} rebuilds every examples/*_example.m.
% FigureScale : double (optional)
%   Name-value argument. Figure pixel density multiplier. Default 2.
%
% Returns
% ---------
% None
%   Writes examples/<name>.mlx and examples/<name>.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    arguments
        names = {}
        opts.FigureScale (1, 1) double {mustBePositive} = 2
    end

    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    examplesFolder = fullfile(repoRoot, 'examples');
    addpath(genpath(fullfile(repoRoot, 'mhkit')));

    if isempty(names)
        files = dir(fullfile(examplesFolder, '*_example.m'));
        names = {files.name};
    end
    names = cellstr(string(names));

    % Without a desktop MATLAB picks a default figure size from the virtual
    % screen, which is much wider than the 560x420 default users see in the
    % desktop. Pin it so figures come out the same size in both.
    set(groot, 'DefaultFigurePosition', [100 100 560 420]);

    % The hidden browser used by the live editor honours Chromium's device
    % scale factor flag
    if opts.FigureScale ~= 1
        scaleFlag = sprintf('--force-device-scale-factor=%g', opts.FigureScale);
        windowManager = matlab.internal.cef.webwindowmanager.instance();
        windowManager.setStartupOptions('InProcess', scaleFlag);
        windowManager.setStartupOptions('ExternalProcess', scaleFlag);
    end

    startDir = pwd;
    cleanup = onCleanup(@() cd(startDir));
    cd(examplesFolder);

    % Force the light theme for this session only (R2025a and newer follow
    % the OS theme by default, which produces dark figures on a dark OS).
    try
        s = settings;
        s.matlab.appearance.MATLABTheme.TemporaryValue = 'Light';
        s.matlab.appearance.figure.GraphicsTheme.TemporaryValue = 'light';
    catch
        % Releases before R2025a have no theme settings and are always light.
    end

    for i = 1:numel(names)
        [~, name] = fileparts(names{i});
        mFile = fullfile(examplesFolder, [name '.m']);
        mlxFile = fullfile(examplesFolder, [name '.mlx']);
        htmlFile = fullfile(examplesFolder, [name '.html']);
        if ~isfile(mFile)
            error('MHKiT:render_example_to_livescript_and_html:NotFound', 'No such example: %s', mFile);
        end

        fprintf('[%d/%d] %s\n', i, numel(names), name);
        fprintf('    converting  %s -> %s\n', [name '.m'], [name '.mlx']);
        matlab.internal.liveeditor.openAndSave(mFile, mlxFile);

        fprintf('    executing   %s\n', [name '.mlx']);
        matlab.internal.liveeditor.executeAndSave(mlxFile);
        close all force;

        fprintf('    exporting   %s (runs the script again)\n', [name '.html']);
        export(mlxFile, htmlFile, 'Run', true);
        close all force;
        setNativeImageSizes(htmlFile, opts.FigureScale);
    end
    fprintf('Done.\n');
end

function setNativeImageSizes(htmlFile, figureScale)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Size HTML figure images to their native size
%
% export() from a hidden browser marks every figure image with
% "width: 100%", which stretches small figures to the page width. This
% replaces that with the figure's native size (pixel size / figureScale)
% so the page looks like a desktop export while keeping the extra pixel
% density.
%
% Parameters
% ------------
% htmlFile : char
%   Path to the exported HTML file, updated in place
% figureScale : double
%   Figure pixel density multiplier used when rendering
%
% Returns
% ---------
% None
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    html = fileread(htmlFile);
    pattern = '<img src="data:image/png;base64,([A-Za-z0-9+/=]+)"[^>]*>';
    [tokens, starts, ends] = regexp(html, pattern, 'tokens', 'start', 'end');
    for k = numel(tokens):-1:1
        png = matlab.net.base64decode(tokens{k}{1}(1:min(64, end)));
        % PNG IHDR: width and height are big-endian uint32 at bytes 17-24
        width = double(typecast(uint8(png(20:-1:17)), 'uint32'));
        height = double(typecast(uint8(png(24:-1:21)), 'uint32'));
        tag = sprintf('<img src="data:image/png;base64,%s" width="%d" height="%d" style="max-width: 100%%; height: auto;">', ...
            tokens{k}{1}, round(width / figureScale), round(height / figureScale));
        html = [html(1:starts(k)-1), tag, html(ends(k)+1:end)];
    end
    fid = fopen(htmlFile, 'w', 'n', 'UTF-8');
    fwrite(fid, html, 'char');
    fclose(fid);
end

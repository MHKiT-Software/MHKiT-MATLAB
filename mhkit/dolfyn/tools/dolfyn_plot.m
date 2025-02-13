function dolfyn_plot(ds, field, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Plot fields from Dolfyn generated data structures
%
% Parameters
% ------------
%     ds: structure
%         Structure from the binary instrument data
%     field: string
%         field from ds to be plotted (ex: 'vel')
%     dim: int or array (optional)
%         for higher dimension fields, which dimension(s) should be plotted
%         can be a single value or array (e.g., 1:3 or [1,3,5])
%     plotArgs: cell array (optional)
%         Additional options to be passed to the plot function
%     cbar_min: double (optional)
%         Minimum value for colorbar in 3D plots
%     cbar_max: double (optional)
%         Maximum value for colorbar in 3D plots
%     width: double (optional)
%         Width of figure in pixels (default: 400 per column)
%     height: double (optional)
%         Height of figure in pixels (default: 300 per row)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Parse input arguments
    arguments
        ds
        field
    end
    arguments (Repeating)
        varargin
    end

    % Parse name-value pairs from varargin
    p = inputParser;
    addParameter(p, 'dim', -1);
    addParameter(p, 'plotArgs', {});
    addParameter(p, 'cbar_min', []);
    addParameter(p, 'cbar_max', []);
    addParameter(p, 'width', []);
    addParameter(p, 'height', []);
    addParameter(p, 'kind', 'plot');  % Default to normal plot
    addParameter(p, 'nbins', 50);     % Default number of bins for histogram
    addParameter(p, 'scale', 'linear');  % Options: 'linear', 'log', 'logx', 'logy'
    parse(p, varargin{:});
    options = p.Results;

    % Basic input validation
    assert(isstruct(ds), ['ds must be a structure generated from one of '...
        'the dolfyn.io functions'])
    assert(isfield(ds,'attrs'), ['It looks like ds was not generated from ' ...
        'one of the dolfyn.io functions'])
    assert(isfield(ds,field), 'The specified field is not present in ds')

    % Get data dimensions
    x_d = length(size(ds.(field).data));
    max_dim = size(ds.(field).data); max_dim = max_dim(end);

    % Handle dimension input
    if x_d > 2
        if isnumeric(options.dim) && length(options.dim) > 1
            dims = options.dim;
        else
            dims = options.dim;
        end
        assert(all(dims >= 1), ['For this field you must specify ' ...
            'which dimension you would like to plot'])
        assert(all(dims <= max_dim), ['The dimension you have ' ...
            'specified is greater than the maximum dimension, %d'], max_dim)
    else
        % For 1D and 2D data, set dims to 1
        dims = 1;
    end

    % Create figure and clear any existing subplots
    if length(dims) > 1
        fig = figure;
        clf(fig);  % Clear the figure
        % Calculate number of rows needed (max 3 plots per row)
        num_plots = length(dims);
        num_cols = min(3, num_plots);
        num_rows = ceil(num_plots / num_cols);

        % Set default width and height if not provided
        if isempty(options.width)
            fig_width = 400 * num_cols;
        else
            fig_width = options.width;
        end

        if isempty(options.height)
            fig_height = 300 * num_rows;
        else
            fig_height = options.height;
        end

        % Create figure with specified size and clear any subplots
        set(fig, 'Position', [100 100 fig_width fig_height]);
    elseif isempty(get(0,'CurrentFigure')) || isempty(get(gcf,'CurrentAxes'))
        % Single plot with specified size
        if ~isempty(options.width) && ~isempty(options.height)
            fig = figure('Position', [100 100 options.width options.height]);
        else
            fig = figure;
        end
        clf(fig);  % Clear the figure
    else
        % If using existing figure, clear it
        clf(gcf);
    end

    % Loop through dimensions if multiple
    for plot_idx = 1:length(dims)
        if length(dims) > 1
            subplot(num_rows, num_cols, plot_idx);
            hold on;  % Enable hold for the subplot
        else
            hold on;  % Enable hold for single plot
        end

        % Main plotting logic
        current_dim = dims(plot_idx);

        if strcmpi(options.kind, 'hist')
            plot_histogram(current_dim);
        else
            switch x_d
                case 2
                    plot_1d_data(current_dim);
                case 3
                    plot_2d_data(current_dim);
                case 4
                    plot_3d_data(current_dim);
            end
        end
    end

    % Adjust subplot spacing if multiple plots
    if length(dims) > 1
        set(gcf, 'Units', 'pixels')
    end

    % Final cleanup to ensure no settings persist
    hold off;  % Release hold state
    if exist('fig', 'var')
        % If we created a new figure, make sure it's clean
        figure(fig);  % Make sure we're on the right figure
        set(gca, 'GridLineStyle', 'none');  % Clear any grid settings
        set(gca, 'XGrid', 'off', 'YGrid', 'off');  % Ensure grid is off
        set(fig, 'NextPlot', 'replace');  % Reset hold state
    else
        % If we used an existing figure, clean it up
        set(gca, 'GridLineStyle', 'none');
        set(gca, 'XGrid', 'off', 'YGrid', 'off');
        hold off;
    end

    % Helper function for histogram plots
    function plot_histogram(current_dim)
        % Get the data based on dimensions
        switch x_d
            case 2
                % 1D data
                plot_data = ds.(field).data;
            case 3
                % 2D data
                plot_data = ds.(field).data(:,:,current_dim);
            case 4
                % 3D data
                plot_data = squeeze(ds.(field).data(:,:,:,current_dim));
        end

        % Flatten the data and remove NaN values
        flat_data = plot_data(:);
        flat_data = flat_data(~isnan(flat_data));

        % Create histogram
        histogram(flat_data, options.nbins);  % Remove normalization

        % Labels
        xlabel(field, 'Interpreter', 'none');
        ylabel('Count', 'Interpreter', 'none');

        % Title
        if x_d > 2
            try
                coord_value = ds.(field).coords.(ds.(field).dims{end}){current_dim};
            catch
                coord_value = string(ds.(field).coords.(ds.(field).dims{end})(current_dim));
            end
            title(strcat(ds.(field).dims{end}, " = ", coord_value));
        else
            if isfield(ds.(field), 'long_name') && ~isempty(ds.(field).long_name)
                title(ds.(field).long_name + " [" + ds.(field).units + "]");
            end
        end

        % Apply scales based on option
        switch options.scale
            case 'log'
                set(gca, 'XScale', 'log', 'YScale', 'log');
            case 'logx'
                set(gca, 'XScale', 'log');
            case 'logy'
                set(gca, 'YScale', 'log');
            case 'linear'
                % Default linear scale, do nothing
            otherwise
                warning('Invalid scale option "%s". Using linear scale.', options.scale);
        end

        % Apply any additional plot arguments
        if ~isempty(options.plotArgs)
            set(gca, options.plotArgs{:});
        end

        hold off;
        set(gca, 'XGrid', 'off', 'YGrid', 'off');
    end

    % Helper function for 1D plots
    function plot_1d_data(current_dim)
        if contains(string(ds.(field).dims{1}), "time")
            x_data = datetime(ds.(field).coords.(ds.(field).dims{1}), ...
                'convertfrom', 'posixtime');
        else
            x_data = ds.(field).coords.(ds.(field).dims{1});
        end
        y_data = ds.(field).data;
        plot(x_data, y_data, options.plotArgs{:})
        xlabel(ds.(field).dims{1}, 'Interpreter', 'none')
        y_label = string(field);
        if isfield(ds.(field), 'units') && ~isempty(ds.(field).units)
            y_label = y_label + " [" + string(ds.(field).units) + "]";
        end
        ylabel(y_label, 'Interpreter', 'none');
        set_title([], current_dim);
        hold off;
        set(gca, 'XGrid', 'off', 'YGrid', 'off');
    end

    % Helper function for 2D plots
    function plot_2d_data(current_dim)
        if contains(string(ds.(field).dims{1}), "time")
            x_data = datetime(ds.(field).coords.(ds.(field).dims{1}), ...
                'convertfrom', 'posixtime');
        else
            x_data = ds.(field).coords.(ds.(field).dims{1});
        end
        y_data = ds.(field).data(:,:,current_dim);
        plot(x_data, y_data, options.plotArgs{:})
        xlabel(ds.(field).dims{1}, 'Interpreter', 'none')
        y_label = string(field);
        if isfield(ds.(field), 'units') && ~isempty(ds.(field).units)
            y_label = y_label + " [" + string(ds.(field).units) + "]";
        end
        ylabel(y_label, 'Interpreter', 'none');
        set_title(ds.(field).coords.(ds.(field).dims{2}){current_dim}, current_dim);
        hold off;
        set(gca, 'XGrid', 'off', 'YGrid', 'off');
    end

    % Helper function for 3D plots
    function plot_3d_data(current_dim)
        if contains(string(ds.(field).dims{1}), "time")
            x_data = datetime(ds.(field).coords.(ds.(field).dims{1}), ...
                'convertfrom', 'posixtime');
        else
            x_data = ds.(field).coords.(ds.(field).dims{1});
        end
        y_dat = ds.(field).coords.(ds.(field).dims{2});
        z_dat = squeeze(ds.(field).data(:,:,:,current_dim))';

        h = pcolor(x_data, y_dat, z_dat);
        shading flat

        % Force axes to be drawn on top
        set(gca, 'Layer', 'top', 'Box', 'on')
        set(gca, 'XGrid', 'off', 'YGrid', 'off')
        set(gca, 'TickDir', 'out')
        set(gca, 'LineWidth', 1)

        xlabel(ds.(field).dims{1}, 'Interpreter', 'none')
        ylabel(strcat(ds.(field).dims{2}, " [m]"), 'Interpreter', 'none')

        % Set colormap based on field type
        switch field
            case 'vel'
                colormap(bluewhitered_colormap(256))
            case 'corr'
                colormap(viridis_colormap(256))
            otherwise
                colormap('default')  % Use MATLAB's default colormap
        end

        % Set colorbar limits if provided
        if ~isempty(options.cbar_min) && ~isempty(options.cbar_max)
            caxis([options.cbar_min options.cbar_max]);
        end

        cb = colorbar;
        cb_label = string(field);
        if isfield(ds.(field), 'units') && ~isempty(ds.(field).units)
            cb_label = cb_label + " [" + string(ds.(field).units) + "]";
        end
        cb.Label.String = cb_label;

        % Set colorbar limits if provided
        if ~isempty(options.cbar_min) && ~isempty(options.cbar_max)
            caxis([options.cbar_min options.cbar_max]);
        end

        try
            coord_value = ds.(field).coords.(ds.(field).dims{3}){current_dim};
        catch
            coord_value = string(ds.(field).coords.(ds.(field).dims{3})(current_dim));
        end
        set_title(coord_value, current_dim);
        hold off;
        set(gca, 'XGrid', 'off', 'YGrid', 'off');
    end

    % Helper function for setting titles
    function set_title(coord_value, current_dim)
        if length(dims) > 1
            if ~isempty(coord_value)
                title(strcat(ds.(field).dims{3}, " = ", coord_value))
            elseif isfield(ds.(field), 'long_name') && ~isempty(ds.(field).long_name)
                title(ds.(field).long_name + " [" + ds.(field).units + "]" + ...
                    " (dim=" + num2str(current_dim) + ")")
            end
        else
            if isfield(ds.(field), 'long_name') && ~isempty(ds.(field).long_name)
                title(ds.(field).long_name + " [" + ds.(field).units + "]")
            end
        end
    end
end

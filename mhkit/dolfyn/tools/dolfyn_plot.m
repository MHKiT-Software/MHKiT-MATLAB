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
%     dim: int (optional)
%         for higher dimension fields, which dimension should be plotted
%     title: string (optional)
%         Title for the plot
%     plotArgs: cell array (optional)
%         Additional options to be passed to the plot function
%
%     call with options -> dolfyn_plot(ds,'vel','dim',1,'title','My plot','plotArgs',{'Color','r'})
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    addParameter(p, 'title', '');
    addParameter(p, 'plotArgs', {});
    parse(p, varargin{:});
    options = p.Results;

    assert(isstruct(ds),['ds must be a structure generated from one of '...
        'the dolfyn.io functions'])

    assert(isfield(ds,'attrs'),['It looks like ds was not generated from ' ...
        'one of the dolfyn.io functions'])

    assert(isfield(ds,field),'The specified field is not present in ds')

    x_d = length(size(ds.(field).data));
    max_dim = size(ds.(field).data); max_dim = max_dim(end);

    if x_d > 2
        assert(options.dim >= 1, ['For this field you must specify ' ...
            'which dimension you would like to plot'])
        assert(options.dim <= max_dim, ['The dimension you have ' ...
            'specified is greater than the maximum dimension, %d'],max_dim)
    end

    if ~isempty(options.title)
        assert(isstring(options.title) || ischar(options.title), ...
            'Title must be character or string')
    end

    switch x_d
        case 2
        % Plot 1D data
        if contains(string(ds.(field).dims{1}),"time")
            x_data = datetime(ds.(field).coords.(ds.(field).dims{1}),...
                'convertfrom','posixtime');
        else
            x_data = ds.(field).coords.(ds.(field).dims{1});
        end
        y_data = ds.(field).data;

        plot(x_data, y_data, options.plotArgs{:})
        xlabel(ds.(field).dims{1},'Interpreter', 'none')
        y_label = string(field);
        if isfield(ds.(field),'units') && ~isempty(ds.(field).units)
            y_label = y_label + " [" + string(ds.(field).units) + "]";
        end
        ylabel(y_label,'Interpreter', 'none');


        case 3
        % Plot 2D data
        if contains(string(ds.(field).dims{1}),"time")
            x_data = datetime(ds.(field).coords.(ds.(field).dims{1}),...
                'convertfrom','posixtime');
        else
            x_data = ds.(field).coords.(ds.(field).dims{1});
        end

        y_data = ds.(field).data(:,:,options.dim);

        plot(x_data, y_data, options.plotArgs{:})
        xlabel(ds.(field).dims{1},'Interpreter', 'none')
        y_label = string(field);
        if isfield(ds.(field),'units') && ~isempty(ds.(field).units)
            y_label = y_label + " [" + string(ds.(field).units) + "]";
        end

        ylabel(y_label,'Interpreter', 'none');

        if ~isempty(options.title)
            title(options.title);
        else
            title(strcat(ds.(field).dims{2}," = ",...
                ds.(field).coords.(ds.(field).dims{2}){options.dim}))
        end

        case 4
        % Plot 3D data
        if contains(string(ds.(field).dims{1}),"time")
            x_data = datetime(ds.(field).coords.(ds.(field).dims{1}),...
                'convertfrom','posixtime');
        else
            x_data = ds.(field).coords.(ds.(field).dims{1});
        end
        y_dat = ds.(field).coords.(ds.(field).dims{2});
        z_dat = squeeze(ds.(field).data(:,:,:,options.dim))';

        surf(x_data, y_dat, z_dat, 'edgecolor','none', options.plotArgs{:})
        view(2)
        xlabel(ds.(field).dims{1},'Interpreter', 'none')
        ylabel(strcat(ds.(field).dims{2}," [m]"),'Interpreter', 'none')
        if ~isempty(options.title)
            title(options.title);
        else
            try
                title(strcat(ds.(field).dims{3}," = ",...
                    ds.(field).coords.(ds.(field).dims{3}){options.dim}))
            catch
                title(strcat(ds.(field).dims{3}," = ", string(...
                    ds.(field).coords.(ds.(field).dims{3})(options.dim))))
            end
        end
        colormap(bluewhitered)
        cb = colorbar;
        cb_label = string(field);
        if isfield(ds.(field),'units') && ~isempty(ds.(field).units)
            cb_label = cb_label + " [" + string(ds.(field).units) + "]";
        end
        cb.Label.String = cb_label;
    end

end

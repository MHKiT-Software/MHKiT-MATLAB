function result = delft_3d_create_points(x, y, waterdepth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Turns three coordinate inputs into a single output struct of points.
%
% Parameters
% ------------
%    x: numeric array or number
%        x values to create points.
%    y: numeric array or number
%        y values to create points.
%    waterdepth: numeric array or number
%        waterdepth values to create points.
%
% Returns
% ---------
%    result: struct
%        A struct containing the points with fields x, y, waterdepth, and df (python output pandas DataFrame).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Input type guards
    if ~isnumeric(x)
        error('MATLAB:delft_3d_create_points:InvalidInput', 'x must be a numeric type.');
    end

    if ~isnumeric(y)
        error('MATLAB:delft_3d_create_points:InvalidInput', 'y must be a numeric type.');
    end

    if ~isnumeric(waterdepth)
        error('MATLAB:delft_3d_create_points:InvalidInput', 'waterdepth must be a numeric type.');
    end

    % Coerce the double into python types
    % The inputs have multiple element convert them into numpy arrays
    if length(x) > 1
        x = py.numpy.array(x);
    end

    if length(y) > 1
        y = py.numpy.array(y);
    end

    if length(waterdepth) > 1
        waterdepth = py.numpy.array(waterdepth);
    end

    % Call Python function to create points
    python_result = py.mhkit.river.io.d3d.create_points(x, y, waterdepth);

    % Convert Python dataframe to struct
    df = python_result;
    result = convert_numeric_dataframe_to_struct(df);
end

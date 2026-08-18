function result = mhkit_restore_spectrum_output(result, input_kind, statistic_name, time)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reassembles a primitive statistic result into the caller's input
% container style
%
% Sister function to mhkit_standardize_spectrum_input. "matrix" and
% "struct" input styles return the primitive numeric result unchanged.
% "table" and "timetable" input styles are reassembled into a table or
% timetable, one row per spectrum, with a single variable named
% statistic_name.
%
% Parameters
% ------------
% result : numeric
%   Primitive statistic result (scalar or column vector, one value per
%   spectrum) from a spectral statistic function
% input_kind : string
%   Value returned by mhkit_standardize_spectrum_input
% statistic_name : string
%   Name for the output variable when input_kind is "table" or
%   "timetable" (e.g. "peak_period")
% time : column vector or [] (optional)
%   Time axis returned by mhkit_standardize_spectrum_input. Required for
%   "timetable"; used as an optional 'time' variable for "table" if
%   present.
%
% Returns
% ---------
% result : numeric, table, or timetable
%   Matches the caller's original input style
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    result {mustBeNumeric}
    input_kind (1,1) string
    statistic_name (1,1) string
    time = []
end

switch input_kind
    case {"matrix", "struct"}
        % unchanged, primitive numeric

    case "table"
        if isempty(time)
            result = table(result(:), 'VariableNames', {char(statistic_name)});
        else
            result = table(time, result(:), 'VariableNames', {'time', char(statistic_name)});
        end

    case "timetable"
        if isempty(time)
            error('MHKiT:mhkit_restore_spectrum_output:InvalidInput', ...
                'time is required to restore a timetable output.');
        end
        result = timetable(time, result(:), 'VariableNames', {char(statistic_name)});

    otherwise
        error('MHKiT:mhkit_restore_spectrum_output:InvalidInput', ...
            'Unrecognized input_kind "%s".', input_kind);
end

end

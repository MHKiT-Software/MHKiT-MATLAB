function names = mhkit_frequency_to_column_names(frequency, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stringifies a frequency vector into valid MATLAB table/timetable
% variable names
%
% Converts each frequency value [Hz] into a self-documenting column name
% (e.g. 0.1 -> "f_0_1000Hz"), for building table/timetable spectrum input
% without resorting to meaningless names like 'f1', 'f2', ...
%
% Parameters
% ------------
% frequency : column vector [Hz]
%   Frequency values, one per table/timetable column
% options.decimals : integer (optional)
%   Number of decimal places included in each name. Default = 4
%
% Returns
% ---------
% names : string array
%   One valid MATLAB identifier per frequency value, e.g. "f_0_1000Hz"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    frequency {mustBeNumeric, mustBeVector, mustBeNonnegative}
    options.decimals (1,1) {mustBeInteger, mustBeNonnegative} = 4
end

arguments (Output)
    names (:,1) string
end

frequency = frequency(:);
formatted = compose("%." + options.decimals + "f", frequency);
names = "f_" + strrep(formatted, ".", "_") + "Hz";

if numel(unique(names)) ~= numel(names)
    error('MHKiT:mhkit_frequency_to_column_names:DuplicateNames', ...
        ['Generated column names are not unique; increase ' ...
        'options.decimals or supply unique frequency values.']);
end

end

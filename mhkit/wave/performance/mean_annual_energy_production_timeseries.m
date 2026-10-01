function maep = mean_annual_energy_production_timeseries(CW, J, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates mean annual energy production (MAEP) from time-series
%
% MAEP = (T/n) * sum(CW * J)
% where T = 8766 hours (average length of a year)
%
% Each sample is weighted equally, so the time series is assumed to be
% representative of a full year of sea states. A warning is issued when
% a time axis is available and the record spans less than one year.
%
% Parameters
% ------------
% CW : vector or timetable [m]
%   Capture width. A timetable must have one variable, and its row
%   times are used as the time axis.
% J : vector or timetable [W/m]
%   Wave energy flux. A timetable must have one variable, and its row
%   times are used as the time axis.
% time : datetime, duration, or numeric vector [s] (optional)
%   Time axis for CW and J, used only to check that the record spans at
%   least one year. Ignored if CW or J is a timetable.
%
% Returns
% ---------
% maep : double [W*h]
%   Mean annual energy production
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    CW {mustBeA(CW, {'numeric', 'timetable'})}
    J {mustBeA(J, {'numeric', 'timetable'})}
    options.time {mustBeA(options.time, {'numeric', 'datetime', 'duration'})} = []
end

arguments (Output)
    maep {mustBeNumeric}
end

[CW_vals, CW_time] = extract_values_and_time(CW, 'CW');
[J_vals, J_time] = extract_values_and_time(J, 'J');

if length(CW_vals) ~= length(J_vals)
    error('MHKiT:mean_annual_energy_production_timeseries:InvalidInput', ...
        'CW length (%d) must match J length (%d)', ...
        length(CW_vals), length(J_vals));
end

% Prefer timetable row times over the optional time argument
if ~isempty(CW_time)
    time = CW_time;
elseif ~isempty(J_time)
    time = J_time;
else
    time = options.time(:);
end

T = 8766;  % Average length of a year in hours per IEC 62600-101 Ed 2.0 Section A.2
n = length(CW_vals);

if ~isempty(time)
    if length(time) ~= n
        error('MHKiT:mean_annual_energy_production_timeseries:InvalidInput', ...
            'time length (%d) must match CW length (%d)', length(time), n);
    end
    span = max(time) - min(time);
    if isnumeric(span)
        span_hours = span / 3600;
    else
        span_hours = hours(span);
    end
    if span_hours < T
        warning('MHKiT:mean_annual_energy_production_timeseries:RecordShorterThanOneYear', ...
            ['Time series spans %.1f days, less than one year (%.2f days). ' ...
             'MAEP weights every sample equally and assumes the record is ' ...
             'representative of a full year.'], span_hours / 24, T / 24);
    end
end

maep = (T / n) * sum(CW_vals(:) .* J_vals(:));

end


function [vals, time] = extract_values_and_time(x, name)
% Returns the data as a column vector plus row times for timetable input.
if istimetable(x)
    if width(x) ~= 1
        error('MHKiT:mean_annual_energy_production_timeseries:InvalidInput', ...
            '%s timetable must have exactly one variable, got %d', name, width(x));
    end
    vals = x{:, 1};
    time = x.Properties.RowTimes;
else
    vals = x;
    time = [];
end
if ~isvector(vals)
    error('MHKiT:mean_annual_energy_production_timeseries:InvalidInput', ...
        '%s must be a vector', name);
end
vals = vals(:);
end

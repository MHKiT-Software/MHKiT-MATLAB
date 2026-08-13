function [spectrum, frequency, time, input_kind, remaining] = mhkit_standardize_spectrum_input(S, function_name, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Standardizes a spectrum given in any supported container to primitive
% numeric arrays
%
% Accepts a PSD spectrum as a plain numeric matrix, a struct-of-vectors,
% a MATLAB table, or a MATLAB timetable, and extracts primitive numeric
% spectrum/frequency/time arrays.
%
% Parameters
% ------------
% S : numeric, struct, table, or timetable
%   numeric: spectrum matrix/vector, frequency required as varargin{1}
%   struct: S.spectrum, S.frequency, optional S.time
%   table: one row per spectrum, one variable (column) per frequency bin;
%     frequency required as varargin{1}; an optional 'time' variable is
%     used as the time axis
%   timetable: one variable per frequency bin (rows = RowTimes), frequency
%     required as varargin{1}
% function_name : string
%   Name of the calling function, used to scope any error raised
% varargin : repeating
%   varargin{1} is frequency when S is numeric, table, or timetable.
%   Any further arguments are returned unconsumed in `remaining`.
%
% Returns
% ---------
% spectrum : matrix
%   Nx K numeric spectrum, frequency as rows, one spectrum per column
% frequency : column vector
%   Nx1 frequency
% time : column vector or []
%   Kx1 datetime/duration if available, otherwise empty
% input_kind : string
%   "matrix", "struct", "table", or "timetable", gets passed to
%   mhkit_restore_spectrum_output
% remaining : cell
%   Unconsumed trailing varargin entries
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    S
    function_name (1,1) string
end

arguments (Repeating)
    varargin
end

time = [];
remaining = varargin;

if istimetable(S)
    input_kind = "timetable";
    if isempty(remaining)
        error(sprintf('MHKiT:%s:InvalidInput', function_name), ...
            '%s requires frequency as an additional argument for timetable input.', function_name);
    end
    frequency = remaining{1};
    remaining(1) = [];
    spectrum = S{:, :}.';
    time = S.Properties.RowTimes;

elseif istable(S)
    input_kind = "table";
    if isempty(remaining)
        error(sprintf('MHKiT:%s:InvalidInput', function_name), ...
            '%s requires frequency as an additional argument for table input.', function_name);
    end
    frequency = remaining{1};
    remaining(1) = [];
    var_names = S.Properties.VariableNames;
    time_idx = strcmpi(var_names, 'time');
    if any(time_idx)
        time = S.(var_names{time_idx});
        S = removevars(S, var_names{time_idx});
    end
    spectrum = S{:, :}.';

elseif isstruct(S)
    input_kind = "struct";
    if ~isfield(S, 'spectrum') || ~isfield(S, 'frequency')
        error(sprintf('MHKiT:%s:InvalidInput', function_name), ...
            '%s requires a struct with spectrum and frequency fields.', function_name);
    end
    spectrum = S.spectrum;
    frequency = S.frequency;
    if isfield(S, 'time')
        time = S.time;
    end

elseif isnumeric(S)
    input_kind = "matrix";
    if isempty(remaining)
        error(sprintf('MHKiT:%s:InvalidInput', function_name), ...
            '%s requires frequency as an additional argument for matrix input.', function_name);
    end
    spectrum = S;
    frequency = remaining{1};
    remaining(1) = [];

else
    error(sprintf('MHKiT:%s:InvalidInput', function_name), ...
        '%s does not support input of class %s.', function_name, class(S));
end

if ~isnumeric(spectrum)
    error(sprintf('MHKiT:%s:InvalidInput', function_name), ...
        '%s requires a numeric spectrum, got %s.', function_name, class(spectrum));
end
if ndims(spectrum) > 2
    error(sprintf('MHKiT:%s:InvalidInput', function_name), ...
        '%s requires a 2-D spectrum (vector or matrix), got %d-D input.', function_name, ndims(spectrum));
end

frequency = mhkit_standardize_user_input_to_column_vectors(frequency, function_name);
mhkit_verify_is_column_vector(frequency, function_name);

if any(frequency < 0)
    error(sprintf('MHKiT:%s:InvalidInput', function_name), ...
        '%s requires non-negative frequency values.', function_name);
end

max_wave_frequency_hz = 100;
if any(frequency > max_wave_frequency_hz)
    warning(sprintf('MHKiT:%s:HighFrequency', function_name), ...
        ['%s received frequency values above %g Hz, which is unusually ' ...
        'high for wave energy spectra - check the frequency units.'], ...
        function_name, max_wave_frequency_hz);
end

% time may be datetime/duration, not numeric, so it can't go through
% mhkit_standardize_user_input_to_column_vectors (numeric-only by design).
if ~isempty(time)
    if isrow(time) && ~isscalar(time)
        time = time.';
    end
    mhkit_verify_is_column_vector(time, function_name);
end

end

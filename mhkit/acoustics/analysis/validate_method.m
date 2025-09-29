function [method_name,method_arg] = validate_method(method)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Helper function for validating the input of possible statistical methods
%
% Parameters
% ------------
%   method: string or cell array
%       Name of statistical method or cell array with method and argument
%
% Returns
% ---------
%   method_name: string
%       Name of validated statistical method
%   method_arg: string or double
%       Corresponding method argument, if applicable
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    method
end

arguments (Output)
    method_name {mustBeText}
    method_arg
end

allowed_methods = [
    "all";
    "any";
    "numel";
    "cumsum";
    "map";
    "max";
    "mean";
    "min";
    "median";
    "prod";
    "quantile";
    "sum";
    "std";
    "var";];

if ~isstring(method) && ~iscell(method)
    error('MHKiT:acoustics:validate_method', "'method' must be a string or a cell array")
end

if isstring(method)
    if ~ismember(method.lower, allowed_methods)
        error('MHKiT:acoustics:validate_method', "'method' is not supported. See list of supported methods.")
    end
    method_name = method.lower;
    method_arg = [];
end

if iscell(method)
    method_name = method{1}.lower;
    if ~ismember(method_name, allowed_methods)
        error('MHKiT:acoustics:validate_method', "'method' is not supported. See list of supported methods")
    end
    method_arg = method{2};
end

end

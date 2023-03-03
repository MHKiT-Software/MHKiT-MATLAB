function out = val_exceeds_thresh(var, options)
% Find values of a variable that exceed a threshold value,
% and assign `options.val` to the velocity data where the threshold is
% exceeded.
%
% Parameters
% ----------
% var : struct (from a fieldname in dataset in create_dataset function)
%   The variable to clean
% options.thresh : double
%   The maximum value of data to screen. Default = 5
% options.val : NaN or double
%   Specifies the value to set the bad values to. Default is `NaN`
%
% Returns
% -------
% out : struct (from a fieldname in dataset in create_dataset function)
%   The input struct `var` with values beyond `options.thresh` set to
%   `options.val`

    % set default arguments
    arguments
        var;
        options.thresh double = 5.0;
        options.val double = NaN;
    end
    
    % copy var.data for output
    out = var;
    % make logical matrix for if elmnt is beyond options.thresh
    flags = abs(out.data) > options.thresh;
    % replace out.data from flags with options.val
    out.data(flags) = options.val;
end


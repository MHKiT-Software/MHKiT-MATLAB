function varargout=magnitude_phase(vector)

%%%%%%%%%%%%%%%%%%%%
%     Calculates magnitude and phase in two or three dimensions
%     call -> [mag, theta] = magnitude_phase({x; y})
%     call -> [mag, theta, phi] = magnitude_phase({x; y; z})
%
% Parameters
% ------------
%     vector: cell array
%           cell array consisting of x, y, and optionally the z component
%           of vector
%       x: array_like x component
%       y: array like y component
%       z: array like z component defined positive up (Optional)
%
%% Returns
% ---------
%     varargout: array
%         depending on number of inputs it will either output an array for
%         magnitude and theta or one for magnitude, theta, and phi
%           magnitude: array - magnitude of vector
%           theta: array - radians from the x axis
%           phi: array - radians from the z-axis defined as positive up
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    vector;
end

% check to see if the input argument is a cell array
if any(~isa(vector, 'cell'))
    ME = MException('MATLAB:magnitude_phase','data must be a cell array ex: {x; y; z}');
    throw(ME);
end

% check how many inputs were given (if not equal to 2 or 3 throw error)
n_output = length(vector);
if n_output ~= 2 && n_output ~= 3
    ME = MException('MATLAB:magnitude_phase','cell array can only be length 2 or 3');
    throw(ME);
end

% make sure number of output variables matches input
if n_output ~= nargout
    ME = MException('MATLAB:magnitude_phase',"Number of output variables doesn't match input format (3 output requested for 2d or 2 output requested for 3d)");
    throw(ME);
end

% Assign cells to arrays
x = vector{1,:};
y = vector{2,:};
three_d = false;
if n_output == 3
    z = vector{3,:};
    three_d = true;
end

if three_d
    if ~isequal(length(x),length(y),length(z))
        ME = MException('MATLAB:magnitude_phase','length of x, y, and z must be equal');
        throw(ME);
    end
    mag = sqrt(x.^2 + y.^2 + z.^2);
    theta = atan2(y,x);
    phi = atan2(sqrt(x.^2 + y.^2),z);
    varargout{1} = mag; varargout{2} = theta; varargout{3} = phi;
else
    if ~isequal(length(x),length(y))
        ME = MException('MATLAB:magnitude_phase','length of x, and y must be equal');
        throw(ME);
    end
    mag = sqrt(x.^2 + y.^2);
    theta = atan2(y,x);
    varargout{1} = mag; varargout{2} = theta;
end


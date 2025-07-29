function cm_data = bluewhitered_colormap(m)
% Creates a diverging colormap from blue to white to red
% The colormap is centered on white and handles the transitions separately
% to avoid hue wrapping issues

% Define the key colors from the reference function
bottom = [0 0 0.5];         % Dark blue
botmiddle = [0 0.5 1];      % Medium blue
middle = [1 1 1];           % White
topmiddle = [1 0 0];        % Medium red
top = [0.5 0 0];           % Dark red

% Create two separate transitions: blue->white and white->red
blue_to_white = [
    bottom;
    botmiddle;
    middle
];

white_to_red = [
    middle;
    topmiddle;
    top
];

if nargin < 1
    % Combine the transitions, removing duplicate white point
    cm_data = [blue_to_white; white_to_red(2:end,:)];
else
    % Calculate half points for interpolation
    m1 = ceil(m/2);
    m2 = floor(m/2);

    % Interpolate each half separately in RGB space
    x1 = linspace(0, 1, size(blue_to_white, 1));
    x2 = linspace(0, 1, size(white_to_red, 1));

    xq1 = linspace(0, 1, m1);
    xq2 = linspace(0, 1, m2);

    % Interpolate directly in RGB space
    first_half = interp1(x1, blue_to_white, xq1);
    second_half = interp1(x2, white_to_red, xq2);

    % Combine the halves, removing duplicate white point
    cm_data = [first_half; second_half(2:end,:)];
end
end

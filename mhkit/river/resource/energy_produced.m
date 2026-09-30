function E = energy_produced(P, seconds)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Returns the energy produced for a given time period provided power
%
% The expected power is computed from a 100-bin histogram of the power
% data, matching MHKiT-Python's use of numpy.histogram and
% scipy.stats.rv_histogram.
%
% Parameters
% ------------
% P : struct
%   Power data
%     P.P : vector or matrix [W]
%       Power
%     P.time : vector [s]
%       Time
% seconds : double [s]
%   Seconds in the time period of interest
%
% Returns
% ---------
% E : double [J]
%   Energy produced in the given length of time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    P struct
    seconds (1,1) {mustBeNumeric}
end

arguments (Output)
    E (1,1) {mustBeNumeric}
end

power_data = P.P(:);

% 100 equal-width bins spanning the data range, matching numpy.histogram,
% which widens a zero-width range by 0.5 on either side
lower = min(power_data);
upper = max(power_data);
if lower == upper
    lower = lower - 0.5;
    upper = upper + 0.5;
end
edges = linspace(lower, upper, 101);
counts = histcounts(power_data, edges);

% Piecewise-constant probability density function of the histogram
density = counts ./ (sum(counts) * diff(edges));

% Evaluate the pdf like scipy.stats.rv_histogram, which assigns x to the
% bin on its right (searchsorted side='right'), so the pdf is 0 at the
% upper edge
x = linspace(edges(1), edges(end), 1000);
bin = discretize(x, edges);
bin(x >= edges(end)) = NaN;
in_range = ~isnan(bin);
pdf_x = zeros(size(x));
pdf_x(in_range) = density(bin(in_range));

% Expected value of power via trapezoidal integration of x*pdf(x)
expected_power = trapz(x, x .* pdf_x);

E = seconds * expected_power;

end

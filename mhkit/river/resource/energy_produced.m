function E=energy_produced(P,seconds)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Returns the energy produced for a given time period provided
%     exceedence probability and power.
%
% Parameters
% ----------
%     P : Power [W]
%
%         structure of form:
%
%           P.P
%
%           P.time [s]
%
%     seconds : float
%         seconds in the time period of interest
%
% Returns
% -------
%     E : float
%         Energy [J] produced in the given length of time
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    P struct
    seconds (1,1) {mustBeNumeric}
end
arguments (Output)
    E (1,1) {mustBeNumeric}
end

power_data = P.P(:);

% Histogram of power with 100 equal-width bins spanning the data range
[counts, edges] = histcounts(power_data, 100);
bin_widths = diff(edges);
total_count = sum(counts);

% Piecewise-constant probability density function of the histogram
density = counts ./ (total_count * bin_widths);

x = linspace(edges(1), edges(end), 1000);
pdf_x = zeros(size(x));

nbins = numel(counts);
for i = 1:nbins
    if i < nbins
        mask = x >= edges(i) & x < edges(i + 1);
    else
        mask = x >= edges(i) & x <= edges(i + 1);
    end
    pdf_x(mask) = density(i);
end

% Expected value of power via trapezoidal integration of x*pdf(x)
expected_power = trapz(x, x .* pdf_x);

E = seconds * expected_power;

end


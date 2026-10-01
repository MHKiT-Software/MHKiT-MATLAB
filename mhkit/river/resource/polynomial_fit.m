function poly = polynomial_fit(x, y, n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Returns a polynomial fit for y given x of order n
%
% Also returns the R-squared score of the fit, computed as the squared
% Pearson correlation between y and the fitted values, matching
% MHKiT-Python's use of scipy.stats.linregress.
%
% Parameters
% ------------
% x : vector
%   x data for polynomial fit
% y : vector
%   y data for polynomial fit, same length as x
% n : integer
%   Order of the polynomial fit
%
% Returns
% ---------
% poly : struct
%   Polynomial fit
%     poly.coef : row vector
%       Polynomial coefficients, highest degree first (as polyval
%       expects)
%     poly.fit : double [-]
%       R-squared coefficient of determination of the fit
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    x {mustBeNumeric, mustBeVector}
    y {mustBeNumeric, mustBeVector}
    n (1,1) {mustBeInteger, mustBeNonnegative}
end

arguments (Output)
    poly struct
end

if numel(x) ~= numel(y)
    error('MHKiT:polynomial_fit:InvalidInput', ...
        'polynomial_fit requires x (%d) and y (%d) to have the same length.', ...
        numel(x), numel(y));
end

% Coefficients ordered highest degree first, matching numpy's poly1d
coef = polyfit(x, y, n);
y_fit = polyval(coef, x);

% R-squared is the squared Pearson correlation between actual and fitted
% values, matching scipy.stats.linregress(y, polynomial_coefficients(x))
correlation = corrcoef(y, y_fit);
r_squared = correlation(1, 2)^2;

poly.coef = coef;
poly.fit = r_squared;

end

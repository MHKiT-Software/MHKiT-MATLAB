function poly=polynomial_fit(x,y,n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Returns a polynomial fit for y given x of order n.
% 
% Parameters
% ----------
%     x : array
%         x data for polynomial fit.
%
%     y : array
%         y data for polynomial fit.
%
%     n : int
%         order of the polynomial fit.
% 
% Returns
% --------
%     poly: structure
%
%
%       poly.coef: polynomial coefficients 
%
%       poly.fit: fit coefficients
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    x {mustBeNumeric, mustBeVector}
    y {mustBeNumeric, mustBeVector}
    n (1,1) {mustBeInteger}
end
arguments (Output)
    poly struct
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



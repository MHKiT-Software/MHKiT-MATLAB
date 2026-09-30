function depth_reg = depth_regime(l, h, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the depth regime based on wavelength and depth
%
% Deep water: h/l > ratio
% This function exists so sinh in wave celerity doesn't blow
% up to infinity.
%
% P.K. Kundu, I.M. Cohen (2000) suggest h/l >> 1 for deep water (pg 209)
% Same citation, they also suggest for 3% accuracy, h/l > 0.28 (pg 210)
% However, since this function allows multiple wavelengths, higher ratio
% numbers are more accurate across varying wavelengths.
%
% Parameters
% ------------
% l : numeric [m]
%   Wavelength (scalar, vector, or array)
% h : double [m]
%   Water column depth
% ratio : double (optional)
%   If h/l > ratio, water depth is classified as deep. Default = 2
%
% Returns
% ---------
% depth_reg : logical
%   Boolean True if deep water, False otherwise.
%   Same shape as input l.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments (Input)
    l {mustBeNumeric}
    h (1,1) {mustBeNumeric, mustBePositive}
    options.ratio (1,1) {mustBeNumeric, mustBePositive} = 2
end

arguments (Output)
    depth_reg logical
end

depth_reg = (h ./ l) > options.ratio;

end

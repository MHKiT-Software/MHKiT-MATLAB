function depth_reg=depth_regime(l,h,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Calculates the depth regime based on wavelength and height
%  Deep water: h/l > ratio
%  This function exists so sinh in wave celerity doesn't blow
%  up to infinity.
%
%  P.K. Kundu, I.M. Cohen (2000) suggest h/l >> 1 for deep water (pg 209)
%  Same citation as above, they also suggest for 3% accuracy, h/l > 0.28 (pg 210)
%  However, since this function allows multiple wavelengths, higher ratio
%  numbers are more accurate across varying wavelengths.
%
% Parameters
% ------------
%    l: vector
%       wave length (m)
%
%    h: double
%         Water depth (m)
%
%    ratio: double or int (optional)
%         If h/l > ratio,
%         water depth will be set to deep. Default ratio = 2
%         to call: energy_flux(k,h,"ratio",ratio)
%
%
% Returns
% -------
%     depth_reg: boolean or boolean array
%         Boolean True if deep water, False otherwise
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    l
    h
    options.ratio = 2;

end

depth_reg=double(py.array.array('d',py.numpy.nditer...
    (py.mhkit.wave.resource.depth_regime(py.numpy.array(l),h,pyargs('ratio',options.ratio)))));


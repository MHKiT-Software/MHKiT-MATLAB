function hs_samples = samples_contour(t_samples,t_contour,hs_contour)
%%%%%%%%%%%%%%%%%%%%%%
%       Get Hs points along a specified environmental contour using
%       user-defined T values.
% 
%       Parameters
%       ----------
%       t_samples : array
%           Points for sampling along return contour
%       t_contour : array
%           T values along contour
%       hs_contour : array
%           Hs values along contour
% 
%       Returns
%       -------
%       hs_samples : array
%           points sampled along return contour

arguments
    t_samples {mustBeNumeric}
    t_contour {mustBeNumeric}
    hs_contour {mustBeNumeric}
end

py.importlib.import_module('mhkit');

t_samples = py.numpy.array(t_samples);
t_contour = py.numpy.array(t_contour);
hs_contour = py.numpy.array(hs_contour);

result = py.mhkit.wave.contours.samples_contour(t_samples, t_contour, hs_contour);

hs_samples = double(result);


end


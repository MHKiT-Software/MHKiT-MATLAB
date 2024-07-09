function environmental_contour=environmental_contours(x1, x2, dt, period, method, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates environmental contours of extreme sea
% states using the improved joint probability distributions
% with the inverse first-order reliability method (IFORM)
% probability for the desired return period (period). Given the
% period of interest a circle of iso-probability is created in the
% in the PCA joint probability (x1, x2) reference frame.
% Using the joint probability value the CDF of the marginal
% distribution is used to find the quantile of each component.
% Finally, using the improved PCA methodology
% the component 2 contour lines are calculated from component 1 using
% the relationships defined in Exkert-Gallup et. al. 2016.
%
% Eckert-Gallup, A. C., Sallaberry, C. J., Dallman, A. R., &
% Neary, V. S. (2016). Application of principal component
% analysis (PCA) and improved joint probability distributions to
% the inverse first-order reliability method (I-FORM) for predicting
% extreme sea states. Ocean Engineering, 112, 307-319.
%
% Parameters
% ------------
%     x1 : vector
%         component 1 data
%
%     x2 : vector
%         component 2 data
%
%     dt : double
%         x1 and x2 sample rate (seconds)
%
%     period : scalar or vector
%         Desired return period (years) for calculation of environmental
%         contour, can be a scalar or a vector.
%
%     PCA: Structure (optional)
% 	      principal component analysis dictionary from previous function
%         call. When supplied the function will skip the PCA calculation
%         for the passe x1, and x2.
%         to call: environmental_contour(x1,x2,dt,period,"PCA",PCA)
%
%     bin_size : double (optional)
%         Data points in each bin
%         to call: environmental_contour(x1,x2,dt,period,"bin_size",bin_size)
%
%     nb_steps : int (optional)
%         Discretization of the circle in the normal space used for
%         IFORM calculation.
%         to call: environmental_contour(x1,x2,dt,period,"nb_steps",nb_steps)
%
%     return_PCA: boolean
% 	      Default False, if True will retun the PCA dictionary
%         to call: environmental_contour(x1,x2,dt,period,"return_PCA",return_PCA)
%
%
% Returns
% ---------
%     environmental_contour: Structure
%         Structure with fields contour1, contour2, and optionally PCA
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments
    x1
    x2
    dt
    period
    method
    options.PCA = py.None;
    options.bin_size = 250;
    options.nb_steps = 1000;
    options.return_fit = py.False;
end

py.importlib.import_module('mhkit');

if options.PCA ~= py.None
   options.PCA = py.dict(options.PCA);
end

if isscalar(period)
    period_py = period;
elseif isvector(period)
    period_py = py.numpy.array(period);
else
    ME = MException('MATLAB:environmental_contour','period must be a vector or scalar');
    throw(ME);
end

data = py.mhkit.wave.contours.environmental_contours(py.numpy.array(x1),py.numpy.array(x2),...
    int32(dt),period_py,method,pyargs('PCA',options.PCA,'bin_size',int32(options.bin_size),...
    'nb_steps',int32(options.nb_steps),'return_fit',options.return_fit));

data_struct = struct(data);

varname1 = strcat(method, "_x1");
varname2 = strcat(method, "_x2");

environmental_contour.contour1 = double(data_struct.(varname1));
environmental_contour.contour2 = double(data_struct.(varname2));

if options.return_fit == py.True
    varfit = strcat(method, "_fit");
    environmental_contour.fit = struct(data_struct.(varfit));
end


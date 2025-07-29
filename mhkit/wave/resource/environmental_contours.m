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
%     method: string or list
%         Copula method to apply. Options include ['PCA','gaussian',
%         'gumbel', 'clayton', 'rosenblatt', 'nonparametric_gaussian',
%         'nonparametric_clayton', 'nonparametric_gumbel', 'bivariate_KDE'
%         'bivariate_KDE_log']
%
%     **OPTIONS
%        min_bin_count: int
%            Passed to _copula_parameters to sets the minimum number of
%            bins allowed. Default = 40.
%        initial_bin_max_val: int, double
%            Passed to _copula_parameters to set the max value of the
%            first bin. Default = 1.
%        bin_val_size: int, double
%            Passed to _copula_parameters to set the size of each bin
%            after the initial bin.  Default 0.25.
%        nb_steps: int
%            Discretization of the circle in the normal space is used for
%            copula component calculation. Default nb_steps=1000.
%        bandwidth:
%            Must specify bandwidth for bivariate KDE method.
%            Default = None.
%        Ndata_bivariate_KDE: int
%            Must specify bivariate KDE method. Defines the contoured
%            space from which samples are taken. Default = 100.
%        max_x1: double
%            Defines the max value of x1 to discretize the KDE space
%        max_x2: double
%            Defines the max value of x2 to discretize the KDE space
%        PCA: dict
%            If provided, the principal component analysis (PCA) on x1,
%            x2 is skipped. The PCA will be the same for a given x1, x2
%            therefore this step may be skipped if multiple calls to
%            environmental contours are made for the same x1, x2 pair.
%            The PCA dict may be obtained by setting return_fit=True when
%            calling the PCA method.
%        return_fit: boolean
%            Will return fitting parameters used for each method passed.
%            Default False.
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
    options.min_bin_count = 40;
    options.initial_bin_max_val = 1;
    options.bin_val_size = 0.25;
    options.nb_steps = 1000;
    options.bandwidth = py.None;
    options.Ndata_bivariate_KDE = 100;
    options.max_x1 = py.None;
    options.max_x2 = py.None;
    options.PCA = py.None;
    options.PCA_bin_size = 250;
    options.return_fit = py.False;
end

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
    int32(dt),period_py,method,pyargs('min_bin_count',int32(options.min_bin_count),...
    'initial_bin_max_val',options.initial_bin_max_val,'bin_val_size',options.bin_val_size,...
    'nb_steps',int32(options.nb_steps),'bandwidth',options.bandwidth,...
    'Ndata_bivariate_KDE',options.Ndata_bivariate_KDE,'max_x1',options.max_x1,...
    'max_x2',options.max_x2,'PCA',options.PCA,'PCA_bin_size',int32(options.PCA_bin_size),...
    'return_fit',options.return_fit));

data_struct = struct(data);

varname1 = strcat(method, "_x1");
varname2 = strcat(method, "_x2");

environmental_contour.contour1 = double(data_struct.(varname1));
environmental_contour.contour2 = double(data_struct.(varname2));

if options.return_fit == py.True
    varfit = strcat(method, "_fit");
    environmental_contour.fit = struct(data_struct.(varfit));
end


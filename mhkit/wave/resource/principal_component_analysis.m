function PCA = principal_component_analysis(x1, x2, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs a modified principal component analysis (PCA)
% [Eckert et. al 2016] on two variables (x1, x2). The additional
% PCA is performed in 5 steps:
% 1) Transform x1 & x2 into the principal component domain and
%    shift the y-axis so that all values are positive and non-zero
% 2) Fit the x1 data in the transformed reference frame with an
%    inverse Gaussian Distribution
% 3) Bin the transformed data into groups of size bin and find the
%    mean of x1, the mean of x2, and the standard deviation of x2
% 4) Perform a first-order linear regression to determine a continuous
%    the function relating the mean of the x1 bins to mean of the
%    x2 bins
% 5) Find a second-order polynomial which best relates the means of
%    x1 to the standard deviation of x2 using constrained
%    optimization
% 
% Eckert-Gallup, A. C., Sallaberry, C. J., Dallman, A. R., &
% Neary, V. S. (2016). Application of principal component
% analysis (PCA) and improved joint probability distributions to
% the inverse first-order reliability method (I-FORM) for predicting
% extreme sea states. Ocean Engineering, 112, 307-319.
% 
% Parameters
% ----------
% x1: array
%     Component 1 data
% x2: array
%     Component 2 data
% bin_size : int
%     Number of data points in each bin
% 
% Returns
% -------
% PCA: structure
%    Fields:
%    -----
%    principal_axes : sign corrected PCA axes
%    shift          : The shift applied to x2
%    x1_fit         : gaussian fit of x1 data
%    mu_param       : fit to _mu_fcn
%    sigma_param    : fit to _sig_fits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments 
    x1
    x2    
    options.bin_size = 250;
end
bin_size = options.bin_size;

assert(isvector(x1),'x1 must be array')
assert(isvector(x2),'x2 must be array')
assert(isfloat(bin_size) || isinteger(bin_size), ...
    'bin_size must be of type int')
bin_size = int32(bin_size);

%% Step 0: Perform Standard PCA
mean_location = 0;
x1_mean_centered = x1 - mean(x1,1);
x2_mean_centered = x2 - mean(x2,1);
n_samples_by_n_features = reshape([x1_mean_centered; x2_mean_centered],...
    length(x1_mean_centered),[]);

[U, ~, Vt] = svd(n_samples_by_n_features,"econ","vector");
Vt = Vt';
[~, max_abs_cols] = max(abs(U));
signs = sign([U(max_abs_cols(1),1),U(max_abs_cols(2),2)]);
principal_axes = Vt .* signs';

%% STEP 1: Transform data into new reference frame
% Apply correct/expected sign convention
principal_axes = abs(principal_axes);
principal_axes(2, 2) = -principal_axes(2, 2);

% Rotate data into Principal direction
x1_and_x2 = reshape([x1; x2],length(x1),[]);
x1_x2_components = x1_and_x2*principal_axes';
x1_components = x1_x2_components(:, 1);
x2_components = x1_x2_components(:, 2);

% Apply shift to Component 2 to make all values positive


end


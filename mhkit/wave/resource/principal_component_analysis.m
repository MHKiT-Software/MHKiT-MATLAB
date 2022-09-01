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
shift = abs(min(x2_components)) + 0.1;
x2_components = x2_components + shift;

%% STEP 2: Fit Component 1 data using a Gaussian Distribution
[x1_sorted, x1_sort_index] = sort(x1_components);
x2_sorted = x2_components(x1_sort_index);

fshape_n = mean(x1_sorted);
x1_fit.scale = length(x1_sorted) / sum(x1_sorted.^(-1) - fshape_n^(-1));
x1_fit.mu = fshape_n / x1_fit.scale;
x1_fit.loc = mean_location;

%% Step 3: Bin Data & find order 1 linear relation between x1 & x2 means
N = length(x1);
minimum_4_bins = floor(N*0.25);
if bin_size > minimum_4_bins
    bin_size = minimum_4_bins;
    fprintf(['To allow for a minimum of 4 bins the bin size has been ' ...
        'set to %d \n'], minimum_4_bins);    
end

N_multiples = floor(N / bin_size);
max_N_multiples_index = N_multiples*bin_size;

x1_integer_multiples_of_bin_size = x1_sorted(1:max_N_multiples_index);
x2_integer_multiples_of_bin_size = x2_sorted(1:max_N_multiples_index);
stride = length(x1_integer_multiples_of_bin_size)/N_multiples;

x1_bins = zeros([stride,N_multiples]);
x2_bins = zeros([stride,N_multiples]);
x1_means  = zeros([N_multiples+1,1]);
x2_means  = zeros([N_multiples+1,1]);
x2_sigmas = zeros([N_multiples+1,1]);

for kk = 1:N_multiples
    jj = (kk-1)*stride + 1;
    qq = (kk)*stride;
    x1_bins(:,kk) = x1_integer_multiples_of_bin_size(jj:qq);
    x2_bins(:,kk) = x2_integer_multiples_of_bin_size(jj:qq);

    x1_means(kk) = mean(x1_bins(:,kk)); 
    x2_means(kk) = mean(x2_bins(:,kk)); 
    x2_sigmas(kk) = std(x2_bins(:,kk));
end

x1_last_bin = x1_sorted(max_N_multiples_index:end);
x2_last_bin = x2_sorted(max_N_multiples_index:end);

x1_means(end) = mean(x1_last_bin); 
x2_means(end) = mean(x2_last_bin); 
x2_sigmas(end) = std(x2_last_bin);

% coeff = polyfit(x1_means, x2_means, 1);
% mu_fit.slope = coeff(1);
% mu_fit.intercept = coeff(2);

mu_fit = linregress(x1_means, x2_means);


%% STEP 4: Find order 2 relation between x1_mean and x2 standard deviation
sigma_polynomial_order = 2;
sig_0 = 0.1 * ones([sigma_polynomial_order + 1,1]);
obj_func = @objective;
con1 = @y_intercept_gt_0;
con2 = @sig_polynomial_min_gt_0;

constraints = struct('con_1',struct('type', 'ineq', 'fun',con1),...
    'con_2',struct('type', 'ineq', 'fun',con2));

sigma_fit = minimize_slsqp( ...
    obj_func, sig_0, x1_means, x2_sigmas, constraints);



    function out = linregress(x,y)
        % Calculate a linear least-squares regression for two sets of 
        % measurements.

        n = length(x);
        xmean = mean(x);
        ymean = mean(y);

        % Average sums of square differences from the mean
        %   ssxm = mean( (x-mean(x))^2 )
        %   ssxym = mean( (x-mean(x)) * (y-mean(y)) )
        dummy = cov(x,y);
        ssxm = dummy(1,1); ssxym = dummy(1,2); ssym = dummy(2,2);

        % R-value
        %   r = ssxym / sqrt( ssxm * ssym )
        if ssxm == 0.0 || ssym == 0.0
            % If the denominator was going to be 0
            r = 0.0;
        else
            r = ssxym / sqrt(ssxm * ssym);
            % Test for numerical error propagation (make sure -1 < r < 1)
            if r > 1.0
                r = 1.0;
            elseif r < -1.0
                r = -1.0;
            end
        end

        slope = ssxym / ssxm;
        intercept = ymean - slope*xmean;

        df = n - 2;  % Number of degrees of freedom
        % n-2 degrees of freedom because 2 has been used up
        % to estimate the mean and standard deviation       
        slope_stderr = sqrt((1 - r^2) * ssym / ssxm / df);

        % Also calculate the standard error of the intercept
        % The following relationship is used:
        %   ssxm = mean( (x-mean(x))^2 )
        %        = ssx - sx*sx
        %        = mean( x^2 ) - mean(x)^2
        intercept_stderr = slope_stderr * sqrt(ssxm + xmean^2);

        out.slope = slope;
        out.intercept = intercept;
        out.rvalue = r;
        out.stderr = slope_stderr;
        out.intercept_stderr = intercept_stderr;
    end

    function error = objective(sig_p, x1_means, x2_sigmas)
        % Polynomial evaluation
        predicted = polyval(sig_p, x1_means);
        % Mean square error
        error = mean((x2_sigmas - predicted).^2);
    end

    function out = y_intercept_gt_0(sig_p)
        out = sig_p(3);
    end

    function out = sig_polynomial_min_gt_0(sig_p)
        out = (sig_p(3) - (sig_p(2)^2) / (4 * sig_p(1)));
    end

end


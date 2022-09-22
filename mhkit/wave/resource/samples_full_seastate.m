function [h_sample, t_sample, weight_points] = samples_full_seastate( ...
    x1, x2, points_per_interval, return_periods, sea_state_duration, ...
    options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample a sea state between contours of specified return periods.
% 
% This function is used for the full sea state approach for the
% extreme load. See Coe et al. 2018 for more details. It was
% originally part of WDRT.
% 
% Coe, R. G., Michelen, C., Eckert-Gallup, A., &
% Sallaberry, C. (2018). Full long-term design response analysis of a
% wave energy converter. Renewable Energy, 116, 356-366.
% 
% Parameters
% ----------
% x1: array
%     Component 1 data
% x2: array
%     Component 2 data
% points_per_interval : int
%     Number of sample points to be calculated per contour interval.
% return_periods: array
%     Vector of return periods that define the contour intervals in
%     which samples will be taken. Values must be greater than zero
%     and must be in increasing order.
% sea_state_duration : int or float
%     x1 and x2 sample rate (seconds)
% method: string or list
%     Copula method to apply. Currently only 'PCA' is implemented.
% bin_size : int
%     Number of data points in each bin
% 
% Returns
% -------
% Hs_Samples: array
%     Vector of Hs values for each sample point.
% Te_Samples: array
%     Vector of Te values for each sample point.
% weight_points: array
%     Vector of probabilistic weights for each sampling point
%     to be used in risk calculations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arguments 
    x1
    x2
    points_per_interval
    return_periods
    sea_state_duration
    options.method = "PCA";
    options.bin_size = 250;
end

method = options.method;
bin_size = options.bin_size;

if ~strcmpi(method,"PCA")
    ME = MException('MATLAB:wave.resource:samples_full_seastate',...
        ['Full sea state sampling is currently only implemented using ' ...
        'the PCA method.']);
    throw(ME);
end
assert(isvector(x1),'x1 must be array')
assert(isvector(x2),'x2 must be array')
assert(isfloat(points_per_interval) || isinteger(points_per_interval), ...
    'points_per_interval must be of type int')
assert(isvector(return_periods),'return_periods must be array')
assert(isfloat(sea_state_duration) || isinteger(sea_state_duration), ...
    'sea_state_duration must be of type int or float')
assert(ischar(method) || isstring(method), 'method must be of type string')
assert(isfloat(bin_size) || isinteger(bin_size), ...
    'bin_size must be of type int')

bin_size = int32(bin_size);

pca_fit = principal_component_analysis(x1, x2, 'bin_size',bin_size);

% Calculate line where Hs = 0 to avoid sampling Hs in negative space
t_zeroline = linspace(2.5, 30, 1000)';
h_zeroline = zeros([length(t_zeroline),1]);

% Transform zero line into principal component space
coeff = pca_fit.principal_axes;
shift = pca_fit.shift;
comp_zeroline = [h_zeroline, t_zeroline]*coeff;
comp_zeroline(:,2) = comp_zeroline(:,2) + shift;

comp1 = pca_fit.x1_fit;
c1_zeroline_prob = ...
    invgauss_cdf(comp_zeroline(:,1), comp1.mu, 0, comp1.scale);

mu_slope = pca_fit.mu_fit.slope;
mu_intercept = pca_fit.mu_fit.intercept;
mu_zeroline = mu_slope * comp_zeroline(:,1) + mu_intercept;

sigma_polynomial_coeffcients = pca_fit.sigma_fit.x;
sigma_zeroline = polyval(sigma_polynomial_coeffcients, comp_zeroline(:,1));
c2_zeroline_prob = cdf(comp_zeroline(:,2),mu_zeroline,sigma_zeroline);

c1_normzeroline = ppf(c1_zeroline_prob, 0, 1);
c2_normzeroline = ppf(c2_zeroline_prob, 0, 1);

contour_probs = 1 ./ ((365*24*60*60)/sea_state_duration .* return_periods);

% Reliability contour generation
% Calculate reliability
beta_lines = ppf((1 - contour_probs), 0, 1);
% Add zero as lower bound to first contour
beta_lines = [0;beta_lines(:)];
% Discretize the circle
theta_lines = linspace(0, 2 * pi, 1000);
% Add probablity of 1 to the reliability set, corresponding to
% probability of the center point of the normal space
contour_probs = [1; contour_probs(:)];

% Vary U1,U2 along circle sqrt(U1^2+U2^2) = beta
u1_lines = transpose(cos(theta_lines).*beta_lines);

% Removing values on the H_s = 0 line that are far from the circles in the
% normal space that will be evaluated to speed up calculations
minval = min(min(u1_lines)) - 0.5;
mask = c1_normzeroline > minval;
c1_normzeroline = c1_normzeroline(mask);
c2_normzeroline = c2_normzeroline(mask);

% Transform to polar coordinates
theta_zeroline = atan2(c2_normzeroline, c1_normzeroline);
rho_zeroline = sqrt(c1_normzeroline.^2 + c2_normzeroline.^2);
theta_zeroline(theta_zeroline < 0) = theta_zeroline(...
        theta_zeroline < 0) + 2 * pi;

[sample_alpha, sample_beta, weight_points] = generate_sample_data(...
        beta_lines, rho_zeroline, theta_zeroline, points_per_interval,...
        contour_probs);


h_sample = 0;
t_sample = 0;
weight_points = 0;


    function out = invgauss_cdf(x, mu, loc, scale)
        % Cumulative distribution function of the given RV.
        log_ndtr = @(x) exp(-x.^2/2);
        x = (x - loc)/scale;
        fac = 1./sqrt(x);
        a = zeros(size(x));
        b = a;
        for i = 1:numel(x)
            a(i) = log(1/sqrt(2*pi)*...
                integral(log_ndtr,-inf,fac(i)*((x(i) / mu) - 1)));
            b(i) = (2/mu) + log(1/sqrt(2*pi)*...
                integral(log_ndtr,-inf,-fac(i)*((x(i) / mu) + 1)));
        end

        out = exp(a + log(1+exp(b-a)));
    end

    function out = cdf(x, loc, scale)
        % Cumulative distribution function of the given RV.
        ndtr = @(x) exp(-x.^2/2);
        out = zeros(size(x));
        x = (x-loc)./scale;
        for i = 1:numel(x)
            out(i) = (1/sqrt(2*pi))*integral(ndtr,-inf,x(i));
        end
    end

    function out = ppf(q, loc, scale)
        % Percent point function (inverse of `cdf`) at q of the given RV.
        out = (-sqrt(2)*erfcinv(2*q));
        out = out .* scale + loc;
    end

    function [sample_alpha, sample_beta, weight_points] = ...
            generate_sample_data(beta_lines, rho_zeroline, ...
            theta_zeroline, points_per_interval, contour_probs)    
        % Calculate radius, angle, and weight for each sample point
        %
        % Parameters
        % ----------
        % beta_lines: np.array
        %     Array of mu fitting function parameters.
        % rho_zeroline: np.array
        %     Array of radii
        % theta_zeroline: np.array
        % points_per_interval: int
        % contour_probs: np.array
        %
        % Returns
        % -------
        % sample_alpha: np.array
        %     Array of fitted sample angle values.
        % sample_beta: np.array
        %     Array of fitted sample radius values.
        % weight_points: np.array
        %     Array of weights for each point.

        num_samples = (length(beta_lines) - 1) * points_per_interval;
        alpha_bounds = zeros([length(beta_lines) - 1, 2]);
        angular_dist = zeros([length(beta_lines) - 1,1]);
        angular_ratio = zeros([length(beta_lines) - 1,1]);
        alpha = zeros([length(beta_lines) - 1, points_per_interval + 1]);
        weight = zeros([length(beta_lines) - 1,1]);
        sample_beta = zeros([num_samples,1]);
        sample_alpha = zeros([num_samples,1]);
        weight_points = zeros([num_samples,1]);

        % Loop over contour intervals
        for i = 1:length(beta_lines)-1
            % Check if any of the radii for the Hs=0, line are smaller than
            % the radii of the contour, meaning that these lines intersect
            r = rho_zeroline - beta_lines(i+1) + 0.01;
            if any(r < 0)
                left = min(find(r < 0));
                right = max(find(r < 0));
                % save sampling bounds
                alpha_bounds(i,:) = [theta_zeroline(left),...
                    theta_zeroline(right) - 2 * pi];
            else
                alpha_bounds(i,:) = [0, 2*pi];
            end
            % Find the angular distance that will be covered by sampling 
            % the disc
            angular_dist(i) = sum(abs(alpha_bounds(i,:)));
            % Calculate ratio of area covered for each contour
            angular_ratio(i) = angular_dist(i) / (2*pi);
            % Discretize the remaining portion of the disc into 10 equally 
            % spaced areas to be sampled
            

        end
    
    end 

end


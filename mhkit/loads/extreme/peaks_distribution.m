classdef peaks_distribution < handle  
    %PEAKS_DISTRIBUTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        qoi_peaks
        method
        paramEsts
        threshold
        npot
        npeaks
    end
    
    methods
        function obj = peaks_distribution()
            obj.method = 'None';
        end

        %% Init methods for ste_block_maxima class
        function weibull(obj, qoi_peaks)
            % Estimate the peaks distribution by fitting a Weibull
            % distribution to the peaks of the response.
            % 
            % The fitted parameters can be accessed through the params
            % field of the returned distribution.
            % 
            % Parameters
            % ----------
            % x : array
            %     Global peaks.
            % 
            % Returns
            % -------
            % parameter estimates set for weibull method.

            obj.qoi_peaks = qoi_peaks;
            obj.method = 'weibull';
            obj.paramEsts = cell([1,4]);

            if ~all(isfinite(obj.qoi_peaks))
                ME = MException('MATLAB:peaks_distribution:weibull',...
                    'The data contains non-finite values.');
                throw(ME);
            end

            obj.paramEsts{1} = 1.0;
            obj.paramEsts{2} = 1.0;
            % Starting point for fit (shape arguments + loc + scale).
            [loc, scale] = obj.fit_loc_scale_support();
            obj.paramEsts{3} = loc;
            obj.paramEsts{4} = scale;
            
            % Use a Nelder-Mead simplex algorithm to find the minimum of 
            % function of one or more variables.
            % 
            % This algorithm has a long history of successful use in 
            % applications. But it will usually be slower than an algorithm 
            % that uses first or second derivative information. In 
            % practice, it can have poor performance in high-dimensional 
            % problems and is not robust to minimizing complicated 
            % functions. Additionally, there currently is no complete 
            % theory describing when the algorithm will successfully
            % converge to the minimum, or how fast it will if it does. 
            % Both the ftol and xtol criteria must be met for convergence.
            vals = obj.minimize_neldermead();
            theta = peaks_distribution.restore([1, 1.0, 0, 1], vals);
            
            if theta(4) < 0 || any(theta(1:2) < 0) 
                ME = MException('MATLAB:peaks_distribution:minimize_neldermead',...
                    ['Optimization converged to parameters that are ' ...
                    'outside the range allowed by the distribution.']);
                throw(ME);
            end

            obj.paramEsts{1} = theta(1); % shape 1 (a)
            obj.paramEsts{2} = theta(2); % shape 2 (c)
            obj.paramEsts{3} = theta(3); % location
            obj.paramEsts{4} = theta(4); % scale
        end

        function weibull_tail_fit(obj, qoi_peaks)
            % Estimate the peaks distribution using the Weibull tail fit
            % method.
            % 
            % The fitted parameters can be accessed through the paramEsts 
            % field ofthe returned distribution.
            % 
            % Parameters
            % ----------
            % x : array
            %     Global peaks.
            % 
            % Returns
            % -------
            % parameter estimates set for weibull tail fit method.

            % Estimate the weibull parameters 
            if isrow(qoi_peaks)
                obj.weibull(qoi_peaks');
            else
                obj.weibull(qoi_peaks);
            end
            obj.method = 'weibull_tail';
            p0 = [obj.paramEsts{2}, obj.paramEsts{4}];

            % Approximate CDF
            x = sort(obj.qoi_peaks);
            npeaks = length(x);
            F = zeros([npeaks,1]);
            for i = 1:npeaks
                F(i) = i / (npeaks + 1.0);
            end

            % Divide into seven sets & fit Weibull
            subset_shape_params = zeros([7,1]);
            subset_scale_params = zeros([7,1]);
            setLim = 0.6:0.05:0.9;
            for qq = 1:7
                xset = x(F > setLim(qq));
                Fset = F(F > setLim(qq));
                [shape, scale, ierr] = obj.minpack_lmdif(xset, Fset, p0);
                if ierr > 3
                    switch ierr                        
                        case 4
                            err_msg = ['the cosine of the angle between ' ...
                                'FVEC and any column of the jacobian is ' ...
                                'at most GTOL in absolute value.'];
                        case 5
                            err_msg = ['number of calls to FCN has reached ' ...
                                'or exceeded MAXFEV.'];
                        case 6
                            err_msg = ['FTOL is too small. No further ' ...
                                'reduction in the sum of squares is ' ...
                                'possible.'];
                        case 7
                            err_msg = ['XTOL is too small. No further ' ...
                                'improvement in the approximate solution ' ...
                                'X is possible.'];
                        case 8
                            err_msg = ['GTOL is too small. FVEC is ' ...
                                'orthogonal to the columns of the jacobian ' ...
                                'to machine precision.'];
                    end
                    ME = MException(['MATLAB:extreme:ste_block_maxima:' ...
                        'minpack_lmdif'],err_msg);
                    throw(ME);
                end
                subset_shape_params(qq) = shape;
                subset_scale_params(qq) = scale;
            end

            obj.paramEsts{1} = 1;                           % shape 1 (a)
            obj.paramEsts{2} = mean(subset_shape_params);   % shape 2 (c)
            obj.paramEsts{3} = 0;                           % location
            obj.paramEsts{4} = mean(subset_scale_params);   % scale
            
        end

        function peaks_over_threshold(obj, qoi_peaks, options)
            % Estimate the peaks distribution using the peaks over 
            % threshold method.
            % 
            % This fits a generalized Pareto distribution to all the peaks 
            % above the specified threshold. The distribution is only 
            % defined for values above the threshold and therefore cannot 
            % be used to obtain integral metrics such as the expected 
            % value. A typical choice of threshold is 1.4 standard 
            % deviations above the mean. The peaks over threshold
            % distribution can be accessed through the `pot` field of the 
            % returned peaks distribution.
            % 
            % Parameters
            % ----------
            % x : array
            %     Global peaks.
            % threshold : float
            %     Threshold value. Only peaks above this value will be used.
            %     Default value calculated as: `np.mean(x) + 1.4 * np.std(x)`
            % 
            % Returns
            % -------
            % parameter estimates set for peaks over threshold method.

            arguments 
                obj
                qoi_peaks
                options.thresh = mean(qoi_peaks) + 1.4*std(qoi_peaks);
            end

            if ~isfloat(options.thresh)
                ME = MException('MATLAB:extreme:peaks_over_threshold',...
                    'threshold must be of type float');
                throw(ME);
            end

            thresh = options.thresh;
            obj.threshold = thresh;
            obj.method = 'pot';
            obj.paramEsts = cell([1,3]);

            % peaks over threshold
            x = sort(qoi_peaks);
            pot = x(x>thresh) - thresh;
            obj.qoi_peaks = pot;
            obj.npot = length(pot);
            obj.npeaks = length(qoi_peaks);

            % Initial guess
            obj.paramEsts{1} = 1.0;
            [loc, scale] = obj.fit_loc_scale_support();
            obj.paramEsts{2} = loc;
            obj.paramEsts{3} = scale;          

            % Fit a generalized Pareto
            vals = obj.minimize_neldermead();
            theta = peaks_distribution.restore([obj.paramEsts{:}], vals);

            if theta(3) < 0 || ~isfinite(theta(1))
                ME = MException('MATLAB:peaks_distribution:minimize_neldermead',...
                    ['Optimization converged to parameters that are ' ...
                    'outside the range allowed by the distribution.']);
                throw(ME);
            end

            obj.paramEsts{1} = theta(1);    % shape
            obj.paramEsts{2} = theta(2);    % location
            obj.paramEsts{3} = theta(3);    % scale   
        end

        %% Statistical Methods
        function out = ppf(obj,q)
            % Percent point function (inverse of cdf) at q of the given 
            % RV.
            % 
            % Parameters
            % ----------
            % q : array_like
            %     lower tail probability
            % 
            % Returns
            % -------
            % out : array_like
            %     quantile corresponding to the lower tail probability q.

            if strcmp(obj.method, 'None')
                ME = MException('MATLAB:extreme:peaks_distribution:ppf',...
                    ['You must initialize the peaks_distribution class with'...
                    'one of the available fit methods.\nMethods include\n' ...
                    '\tweibull\n\tweibull_tail_fit\n\tpeaks_over_threshold'...
                    ' (ppf not available for this method']);
                throw(ME);
            end

            if strcmp(obj.method, 'pot')
                ME = MException('MATLAB:extreme:peaks_distribution:ppf',...
                    ['Percent point function is not available for the ' ...
                    'peaks over threshold method']);
                throw(ME);
            end
            
            [args, loc, scale] = obj.als();

            cond0 = scale > 0;
            cond1 = (0 < q) & (q < 1);
            cond = cond0 & cond1;

            if any(cond)
                a = args(1); c = args(2);                
                out = (-1*log(-q.^(1.0/a)+1)).^(1.0/c);
                out = out * scale + loc;
            end

            if isempty(out)
                out = [];
            end
        end

        function out =  pdf(obj, x)
            % Probability density function at x of the given RV.
            % 
            % Parameters
            % ----------
            % x : array_like
            %     quantiles
            % 
            % Returns
            % -------
            % pdf : array
            %     Probability density function evaluated at x

            if strcmp(obj.method, 'None')
                ME = MException('MATLAB:extreme:peaks_distribution:pdf',...
                    ['You must initialize the peaks_distribution class with'...
                    'one of the available fit methods.\nMethods include\n' ...
                    '\tweibull\n\tweibull_tail_fit\n\t' ...
                    'peaks_over_threshold']);
                throw(ME);
            end

            if strcmp(obj.method, 'pot')
                args = [];
                loc = 0;
                scale = 1;
            else
                [args, loc, scale] = obj.als();
            end

            x = (x-loc)/scale;
            cond0 = scale > 0;
            cond1 = 0 <= x & x <= inf & scale > 0;
            cond = cond0 & cond1;
            out = zeros(size(cond));

            if any(cond)
                goodargs = x(cond);
                if length(args) < 2
                    % Find the nth derivative of a function at a point.
                    % Given a function, use a central difference formula 
                    % with spacing dx to compute the nth derivative at x0.
                    order = 5;
                    dx = 1e-5;
                    weights = [1,-8,0,8,-1]/12.0;
                    val = 0;
                    ho = bitshift(order,-1);
                    for k = 1:order
                        val = val + weights(k)*obj.cdf(x+(k-ho)*dx);
                    end
                    out = val./dx;    
                else
                    negxc = -goodargs.^args(2);
                    exm1c = -1*(exp(negxc)-1);
                    logpdf = (log(args(1)) + log(args(2)) + ...
                        (args(1)-1.0)*log(exm1c)...
                        + negxc + (args(2)-1.0)*log(goodargs));
                    out(cond) = exp(logpdf)./scale;
                end
            end              
        end

        function out =  cdf(obj, x, arg)
            % Cumulative distribution function of the given RV.
            % 
            % Parameters
            % ----------
            % x : array_like
            %     quantiles
            %
            % Returns
            %    -------
            % cdf : ndarray
            %       Cumulative distribution function evaluated at x

            if strcmp(obj.method, 'None')
                ME = MException('MATLAB:extreme:peaks_distribution:cdf',...
                    ['You must initialize the peaks_distribution class with'...
                    'one of the available fit methods.\nMethods include\n' ...
                    '\tweibull\n\tweibull_tail_fit\n\t' ...
                    'peaks_over_threshold']);
                throw(ME);
            end

            if strcmp(obj.method, 'pot')
                out = zeros(size(x));
                out(x < obj.threshold) = nan;
                xt = x(x > obj.threshold);               
                    
                if ~isempty(xt)
                    cdf_inp = xt - obj.threshold;
                    %
                    [args, loc, scale] = obj.als();
                    cdf_inp = (cdf_inp-loc)/scale;
                    cond0 = scale > 0;                    
                    cond1 = peaks_distribution.support_mask(...
                            cdf_inp,args) & scale > 0;                        
                    cond = cond0 & cond1;
                    pot_ccdf = ones(size(cond));        
                    if any(cond)
                        goodargs = -cdf_inp(cond);
                        c = -args(1); 
                        if c == 0
                            bxcx = -1*(exp(goodargs)-1);
                        else
                            bxcx = ...
                                -1*(exp((log(c.*goodargs + 1)./c))-1);
                        end
                        pot_ccdf(cond) = bxcx;
                    end 
                    %
                    pot_ccdf = 1 - pot_ccdf;
                    prop_pot = obj.npot/obj.npeaks;
                    out(x >= obj.threshold) = ...
                        1.0 - (prop_pot * pot_ccdf);
                end
                
            else
                if nargin < 3 
                    [args, loc, scale] = obj.als();
                else
                    loc = arg(3);
                    scale = arg(4);
                    args = arg(1:2);
                end
    
                x = (x-loc)/scale;
                cond0 = scale > 0;
                cond1 = (0 < x) & (x < inf) & scale > 0;
                cond = cond0 & cond1;
                out = zeros(size(cond));
    
                if any(cond)
                    goodargs = x(cond);
                    a = args(1); c = args(2); 
                    exm1c = -1*(exp(-goodargs.^c)-1);
                    out(cond) = exm1c.^a;
                end 
            end             
        end

        function out = expect(obj)
            % Calculate expected value of a function with respect to the 
            % distribution for discrete distribution by numerical summation

            % The expected value of a function ``f(x)`` with respect to a
            % distribution ``dist`` is defined as::
            % 
            %             ub
            %     E[f(x)] = Integral(f(x) * dist.pdf(x)),
            %             lb
            % 
            % where ub and lb are arguments and x has the dist.pdf(x)
            % distribution. If the bounds lb and ub correspond to the
            % support of the distribution, e.g. [-inf, inf] in the default
            % case, then the integral is the unrestricted expectation of 
            % f(x). Also, the function f(x) may be defined such that f(x) 
            % is 0 outside a finite interval in which case the expectation 
            % is  calculated within the finite range [lb, ub].

            if strcmp(obj.method, 'None')
                ME = MException('MATLAB:extreme:peaks_distribution:excpect',...
                    ['You must initialize the peaks_distribution class with'...
                    'one of the available fit methods.\nMethods include\n' ...
                    '\tweibull\n\tweibull_tail_fit\n\tpeaks_over_threshold'...
                    ' (expect not available for this method']);
                throw(ME);
            end

            if strcmp(obj.method, 'pot')
                ME = MException('MATLAB:extreme:peaks_distribution:expect',...
                    ['Expect function is not available for the ' ...
                    'peaks over threshold method']);
                throw(ME);
            end

            [args, loc, scale] = obj.als();
            
            a = 0.0;
            b = inf;
            lb = loc + a * scale;
            ub = loc + b * scale;

            out = obj.quad(lb,ub);
        end   
    
    end

    methods (Access='private')
        %% Fit methods
        function [loc_hat, scale_hat] = fit_loc_scale_support(obj)
            % Estimate loc and scale parameters from data accounting 
            % for support.
            %
            % Parameters
            % ----------           
            %
            % Returns
            % -------
            % loc_hat : float
            %     Estimated location parameter for the data.
            % scale_hat : float
            %     Estimated scale parameter for the data.
            %
            % Estimate location and scale according to the method of 
            % moments.
            [loc_hat, scale_hat] = obj.fit_loc_scale();

            % Compute the support according to the shape parameters.            
            a = 0.0;
            b = inf;
            support_width = b - a;

            %If the support is empty then return the moment-based estimates
            if support_width <= 0
                return
            end

            % Compute the proposed support according to the loc and scale
            % estimates.
            a_hat = loc_hat + a * scale_hat;
            b_hat = loc_hat + b * scale_hat;

            % Use the moment-based estimates if they are compatible with 
            % the data.
            data_a = min(obj.qoi_peaks);
            data_b = max(obj.qoi_peaks);
            if a_hat < data_a && data_b < b_hat
                return 
            end

            % Otherwise find other estimates that are compatible with 
            % the data.
            data_width = data_b - data_a;
            rel_margin = 0.1;
            margin = data_width * rel_margin;

            % For a finite interval, both the location and scale
            % should have interesting values.
            if support_width < inf
                loc_hat = (data_a - a) - margin;
                scale_hat = (data_width + 2 * margin) / support_width;
                return
            end

            % For a one-sided interval, use only an interesting 
            % location parameter
            if a > -inf
                loc_hat = (data_a - a) - margin;
                scale_hat = 1;
                return
            elseif b < inf
                loc_hat = (data_b - b) + margin;
                scale_hat = 1;
                return
            else
                ME = MException('MATLAB:extreme:ste_block_maxima:fit_loc_scale_support',...
                    'Runtime Error');
                throw(ME);
            end
        end

        function [loc_hat, scale_hat] = fit_loc_scale(obj)
            % Estimate loc and scale parameters from data accounting 
            % for support.
 
            % Parameters
            % ----------
            % data : array_like
            %     Data to fit.
     
            % Returns
            % -------
            % loc_hat : float
            %     Estimated location parameter for the data.
            % scale_hat : float
            %     Estimated scale parameter for the data.
            [mu, mu2] = obj.stats();

            muhat = mean(obj.qoi_peaks);
            % Variance (matlab's var() function does not return the correct
	        % value
            mu2hat = mean(abs(obj.qoi_peaks -...
		        mean(obj.qoi_peaks)).^2);

            scale_hat = sqrt(mu2hat / mu2);
            loc_hat = muhat - scale_hat*mu ;

            if isinf(loc_hat) || isnan(loc_hat)
                loc_hat = 0;
            end
            if ~(~isinf(scale_hat) && 0 < scale_hat) || isnan(scale_hat)
                scale_hat = 1;
            end
        end

        function min_params = minimize_neldermead(obj)
            % Minimization of scalar function of one or more variables using the
            % Nelder-Mead algorithm.
            %             
            % References
            % ----------
            % .. [1] Gao, F. and Han, L.
            %    Implementing the Nelder-Mead simplex algorithm with adaptive
            %    parameters. 2012. Computational Optimization and Applications.
            %    51:1, pp. 259-277

            xatol = 1e-4;
            fatol = 1e-4;
            rho = 1;
            chi = 2;
            psi = 0.5;
            sigma = 0.5;
            nonzdelt = 0.05;
            zdelt = 0.00025;

            x0 = [1.0, 1.0];   
            args = obj.als();
            N = length(x0);
            sim = zeros([N + 1, N]);
            sim(1,:) = x0;
            for k = 1:N
                y = x0;
                if y(k) ~= 0
                    y(k) = (1+ nonzdelt)*y(k);
                else
                    y(k) = zdelt;
                end
                sim(k+1,:) = y;
            end

            maxiter = N * 200;

            one2np1 = 2:N+1;
            fsim = zeros([N + 1,1]);

            for k = 1:N+1                
                theta = peaks_distribution.restore(...
                    [args, 0, 1], sim(k,:));
                fsim(k,1) = obj.penalize_nnlf(theta);
            end

            [fsim, ind] = sort(fsim);
            sim = sim(ind,:);

            iterations = 1;

            while iterations < maxiter
                if max(reshape(abs((sim(2:end,:) - sim(1,:))).',1,[])) ...
                        <= xatol && max(abs(fsim(1)-fsim(2:end))) <= fatol
                    break;
                end

                xbar = sum(sim(1:end-1,:),1)/N;
                xr = (1 + rho) * xbar - rho * sim(end,:);
                theta = ...
                    peaks_distribution.restore([args, 0, 1], xr);
                fxr = obj.penalize_nnlf(theta);
                doshrink = false;

                if fxr < fsim(1)
                    xe = (1 + rho * chi) * xbar - rho * chi * sim(end,:);
                    theta = ...
                    peaks_distribution.restore([args, 0, 1], xe);
                    fxe = obj.penalize_nnlf(theta);

                    if fxe < fxr
                        sim(end,:) = xe;
                        fsim(end) = fxe;
                    else
                        sim(end,:) = xr;
                        fsim(end) = fxe;
                    end
                else % fsim(1) <= fxr
                    if fxr < fsim(end-1)
                        sim(end,:) = xr;
                        fsim(end) = fxr;
                    else %fxr >= fsim(end-1)
                        % Perform contraction
                        if fxr < fsim(end)
                            xc = (1 + psi*rho)*xbar - psi*rho*sim(end,:);
                            theta = peaks_distribution.restore(...
                                [args, 0, 1], xc);
                            fxc = obj.penalize_nnlf(theta);

                            if fxc <= fxr
                                sim(end,:) = xc;
                                fsim(end) = fxc;
                            else
                                doshrink = true;
                            end
                        else
                            % Perform an inside contraction
                            xcc = (1 - psi) * xbar + psi * sim(end,:);
                            theta = ...
                            peaks_distribution.restore(...
                            [args, 0, 1], xcc);
                            fxcc = obj.penalize_nnlf(theta);

                            if fxcc < fsim(end)
                                sim(end,:) = xcc;
                                fsim(end) = fxcc;
                            else
                                doshrink = true;
                            end
                        end

                        if doshrink
                            for j = one2np1
                                sim(j,:) = sim(1,:) + sigma*(sim(j,:) - ...
                                    sim(1,:));
                                theta = peaks_distribution.restore(...
                                    [args, 0, 1], sim(j,:));
                                fsim(j,:) = obj.penalize_nnlf(theta);
                            end
                        end
                    end                    
                end
                [fsim, ind] = sort(fsim);
                sim = sim(ind,:);
                iterations = iterations + 1;
            end % end while

            min_params = sim(1,:);
        end

        function out = penalize_nnlf(obj, theta)
            % Penalized negative loglikelihood function.
            % 
            % i.e., - sum (log pdf(x, theta), axis=0) + penalty
            % where theta are the parameters (including loc and scale)
            loc = theta(end-1);
            scale = theta(end);
            args = theta(1:end-2);
            if scale <= 0
                out = inf;
                return
            end
            x = (obj.qoi_peaks - loc)./scale;
            n_log_scale = length(x) * log(scale);

            cond0 = ~((0 < x) & (x < inf));
            n_bad = nnz(cond0);
            if n_bad > 0
                x = x(~cond0);
            end
            if length(args) < 2
                if args == 0
                    logpdf = -x;
                else
                    logpdf = -1*((args+1)*log(1+ args.*x))./args;
                    logpdf(imag(logpdf)~=0) = nan;
                end
            else
                negxc = -x.^args(2);
                exm1c = -1*(exp(negxc)-1);
                logpdf = (log(args(1)) + log(args(2)) + ...
                    (args(1)-1.0)*log(exm1c)...
                    + negxc + (args(2)-1.0)*log(x));
            end
            finite_logpdf = isfinite(logpdf);
            n_bad = n_bad + nnz(~finite_logpdf);
            if n_bad > 0
                pentalty = n_bad * log(1.7976931348623157e+308)*100.0;
                out = -1*sum(logpdf(finite_logpdf)) + pentalty +...
                    n_log_scale;
            else
                out = -1*sum(logpdf) + n_log_scale;
            end
        end

        function [shape, scale, info] = minpack_lmdif(obj, x_data, y_data, x)
            % LMDIF minimizes M functions in N variables by the 
            % Levenberg-Marquardt method.
            %
            %  Discussion:
            %
            %    LMDIF minimizes the sum of the squares of M nonlinear 
            %    functions in N variables by a modification of the 
            %    Levenberg-Marquardt algorithm. The user must provide a 
            %    subroutine which calculates the functions. The jacobian is 
            %    then calculated by a forward-difference approximation.
            %
            %  Licensing:
            %
            %    This code may freely be copied, modified, and used for any 
            %    purpose.           
            %
            %  Author:
            %
            %    Original FORTRAN77 version by Jorge More, Burton Garbow, 
            %    Kenneth Hillstrom. FORTRAN90 version by John Burkardt.
            %
            %  Reference:
            %
            %    Jorge More, Burton Garbow, Kenneth Hillstrom,
            %    User Guide for MINPACK-1,
            %    Technical Report ANL-80-74,
            %    Argonne National Laboratory, 1980.
            %
            %  Parameters:
            %
            %    Input, external FCN, the name of the user-supplied 
            %    subroutine which calculates the functions. The routine 
            %    should have the form:
            %      subroutine fcn ( m, n, x, fvec, iflag )
            %      integer m
            %      integer n
            %      real ( kind = rk ) fvec(m)
            %      integer iflag
            %      real ( kind = rk ) x(n)            
            %
            %    Input/output DIAG(N).  If MODE = 1, then DIAG is set
            %    internally.  If MODE = 2, then DIAG must contain positive 
            %    entries that serve as multiplicative scale factors for the 
            %    variables.
            %
            %    Input, integer MODE, scaling option.
            %    1, variables will be scaled internally.
            %    2, scaling is specified by the input DIAG vector.
            %
            %    Input, FACTOR, determines the initial step bound.  
            %    This bound is set to the product of FACTOR and the 
            %    euclidean norm of DIAG*X if nonzero, or else to FACTOR 
            %    itself.  In most cases, FACTOR should lie in the interval 
            %    (0.1, 100) with 100 the recommended value.
            %
            ftol = 1.49012e-08;
            xtol = 1.49012e-08;
            gtol = 0.0;
            maxfev = 5000;
            max_check = -1;
            n = length(x);
            mode = 1;
            info = 0;
            diag = zeros([1,n]);
            qtf = zeros([1,n]);
            factor = 100;

            fun = @(x,y,args) obj.cdf(x, args)-y;
            % Evaluate the function at the starting point and calculate 
            % its norm.
            args = [1, x(1), 0, x(2)];
            fvec = fun(x_data, y_data, args);
            nfev = 1;
            fnorm = norm(fvec);

            % Initialize Levenberg-Marquardt parameter and iteration 
            % counter.
            par = 0.0;
            iter = 1;

            % Beginning of the outer loop.
            while true
                % Calculate the jacobian matrix.
                fjac = peaks_distribution.fdjac2(...
                    fun, x, fvec, x_data, y_data);
                nfev = nfev + 1;

                % Compute the QR factorization of the jacobian.
                pivot = true;
                [fjac, ipvt, wa1, wa2] = ...
                    peaks_distribution.qrfac(fjac, pivot);

                % On the first iteration and if MODE is 1, scale according
                % to the norms of the columns of the initial jacobian.
                if iter == 1
                    if mode ~= 2
                        diag(1,1:n) = wa2(1,1:n);
                        for j = 1:n
                            if wa2(1,j) == 0.0
                                diag(1,j) = 1.0;
                            end
                        end
                    end

                    % On the first iteration, calculate the norm of the 
                    % scaled X and initialize the step bound DELTA.
                    wa3 = diag .* x;
                    xnorm = sqrt(sum(wa3.^2));
                    delta = factor * xnorm;
                    if delta == 0
                        delta = factor;
                    end                    
                end

                % Form Q' * FVEC and store the first N components in QTF
                wa4 = fvec;
                for j = 1:n
                    if fjac(j,j) ~= 0.0
                        sum2 = dot(wa4(j:end),fjac(j:end,j));
                        temp = -sum2/fjac(j,j);
                        wa4(j:end) = wa4(j:end) + fjac(j:end,j)*temp;
                    end
                    fjac(j,j) = wa1(j);
                    qtf(j) = wa4(j);
                end

                % Compute the norm of the scaled gradient 
                gnorm = 0.0;
                if fnorm ~= 0.0
                    for j = 1:n
                        l = ipvt(j);
                        if wa2(l) ~= 0.0
                            sum2 = 0.0;
                            for i = 1:j
                                sum2 = sum2 + fjac(i,j) * (qtf(i)/fnorm);
                            end
                            gnorm = max(gnorm,abs(sum2/wa2(l)));
                        end
                    end
                end

                %Test for convergence of the gradient norm.
                if gnorm < gtol
                    info = 4;
                    break;
                end

                % Rescale if necessary.
                if mode ~= 2
                    for j = 1:n
                        diag(j) = max(diag(j), wa2(j));
                    end
                end

                % Beginning of inner loop
                while true
                    % Determine the Levenberg-Marquardt parameter.
                    [fjac, par, wa1,  wa2] = peaks_distribution.lmpar(...
                        fjac, ipvt, diag, qtf, delta, par);

                    % Store the direction P and X + P. 
                    % Calculate the norm of P.
                    wa1(1:n) = -wa1(1:n);
                    wa2(1:n) = x(1:n) + wa1(1:n);
                    wa3(1:n) = diag(1:n) .* wa1(1:n);

                    pnorm = norm(wa3);

                    % On the first iteration, adjust the initial step bound
                    if iter == 1
                        delta = min(delta, pnorm);
                    end

                    % Evaluate the function at X + P and calculate its norm
                    args = [1, wa2(1), 0, wa2(2)];
                    wa4 = fun(x_data, y_data, args);
                    nfev = nfev + 1;
                    fnorm1 = norm(wa4);

                    % Compute the scaled actual reduction.
                    if 0.1*fnorm1 < fnorm
                        actred = 1.0 - (fnorm1/fnorm)^2;
                    else
                        actred = -1.0;
                    end

                    % Compute the scaled predicted reduction and the scaled 
                    % directional derivative.
                    for j = 1:n
                        wa3(j) = 0.0;
                        l = ipvt(j);
                        temp = wa1(l);
                        wa3(1:j) = wa3(1:j) + fjac(1:j,j)' * temp;                        
                    end

                    temp1 = norm(wa3)/fnorm;
                    temp2 = (sqrt(par)*pnorm)/fnorm;
                    prered = temp1^2 + temp2^2 / 0.5;
                    dirder = -1*(temp1^2 + temp2^2);

                    % Compute the ratio of the actual to the predicted 
                    % reduction.
                    ratio = 0.0;
                    if prered ~= 0.0
                        ratio = actred / prered;
                    end

                    % Update the step bound.
                    if ratio <= 0.25
                        if actred >= 0.0
                            temp = 0.5;
                        end

                        if actred < 0.0
                            temp = 0.5 * dirder/(dirder + 0.5*actred);
                        end

                        if 0.1*fnorm1 >= fnorm || temp < 0.1
                            temp = 0.1;
                        end
                    else
                        if par == 0.0 || ratio >= 0.75
                            delta = 2.0 * pnorm;
                            par = 0.5*par;
                        end
                    end

                    % Test for successful iteration.

                    if 0.0001 <= ratio
                        x(1:n) = wa2(1:n);
                        wa2(1:n) = diag(1:n) .* x(1:n);
                        fvec = wa4;
                        xnorm = norm(wa2);
                        fnorm = fnorm1;
                        iter = iter + 1;
                    end

                    % Tests for convergence.
                    if abs(actred) <= ftol && prered <= ftol && ...
                            0.5*ratio <= 1.0 
                        info = 1;
                    end

                    if delta <= xtol*xnorm
                        info = 2;
                    end

                    if abs(actred) <= ftol && prered <= ftol && ...
                            0.5*ratio <= 1.0 && info == 2
                        info = 3;                        
                    end

                    if info ~= 0
                        shape = x(1);
                        scale = x(2);
                        return;
                    end

                    % Tests for termination and stringent tolerances.
                    if nfev >= maxfev
                        % it may just be converging slowly so check if the
                        % error is progressing downward and up maxfev if it
                        % seems to be making progress
                        if max_check < 0
                            % first time hitting this mark
                            if ftol/actred > 0.5
                                maxfev = maxfev * 1.5;
                                max_check = actred;
                            else
                                info = 5;
                            end
                        else
                            % its been increased already. check to make
                            % sure its still making progress
                            if actred/max_check < 0.5
                                maxfev = maxfev + 500;
                                max_check = actred;
                            else
                                info = 5;
                            end
                        end
                    end

                    if abs(actred) <= eps && prered <= eps &&...
                            0.5*ratio <= 1.0
                        info = 6;
                    end

                    if delta <= eps*xnorm
                        info = 7;
                    end

                    if gnorm <= eps
                        info = 8;
                    end

                    if info ~= 0
                        shape = x(1);
                        scale = x(2);
                        return;
                    end

                    % End of the inner loop.  
                    % Repeat if iteration unsuccessful.
                    if ratio > 0.0001
                        break; % go to outer loop
                    end

                end % End inner loop
            end % End outer loop


        end

        %% Utility Functions
        function [mu, mu2] = stats(obj)
            % statistics of the given RV using mean variance.            
            [pos, ~,~] = obj.als;
            loc = 0;
            scale = 1; 
            if length(pos) < 2
                obj.paramEsts{2} = loc;
                obj.paramEsts{3} = scale;
                mu = inf;
                mu2 = nan;
            else
                obj.paramEsts{3} = loc;
                obj.paramEsts{4} = scale;
                
                % mu is integral of ppf between 0 and 1
                fun = @(x) obj.ppf(x); 
                mu = integral(fun,0,1);
                fun = @(x) obj.ppf(x).^2;
                mu2p = integral(fun,0,1);
                mu = mu * scale + loc;
                mu2 = mu2p - mu^2;
                mu2 = mu2 * scale * scale;
            end         
        end

        function [args, loc, scale] = als(obj)
            if contains(obj.method,'weibull')
                args = [obj.paramEsts{1},obj.paramEsts{2}];
                loc = obj.paramEsts{3};
                scale = obj.paramEsts{4};
            else
                args = obj.paramEsts{1};
                loc = obj.paramEsts{2};
                scale = obj.paramEsts{3};
            end
        end
                
        function out = quad(obj, a, b)
            % Compute a definite integral.
            % 
            % Integrate func from `a` to `b` (possibly infinite interval) 
            % using a technique from the Fortran library QUADPACK.
            
            % check the limits of integration: \int_a^b, expect a < b
            flip = b < a;
            a = min(a, b);
            b = max(a, b);

            if (a==Inf && b==Inf) || (a==-Inf && b==-Inf)
                ME = MException('MATLAB:extreme:ste_block_maxima:quad',"Infinity " + ...
                    "comparisons don't work with this method.");
                throw(ME);
            end
            
            fun = @(x) x .*obj.pdf(x);
            out = integral(fun,a,b);

            if flip
                out = -out;
            end            
        end
        
    end

    methods(Static, Access='private')
        function new_theta = restore(args, theta)
            i = 1;
            if length(args) < 4
                fixedn = [2];
            else
                fixedn = [1, 3];
            end
            for qq = 1:length(args)
                if ~ismember(qq,fixedn)
                    args(qq) = theta(i);
                    i = i + 1;
                end
            end
            new_theta = args;
        end
        
        function fjac = fdjac2(fcn, x, fvec, x_data, y_data)
            % fdjac2 estimates an M by N jacobian matrix using forward 
            % differences.
            %
            %  Discussion:
            %
            %    This function computes a forward-difference approximation
            %    to the M by N jacobian matrix associated with a specified
            %    problem of M functions in N variables.
            %
            %  Parameters:
            %
            %    FCN, the name of the user-supplied subroutine which
            %    calculates the functions.  The routine should have the 
            %    form:
            %      subroutine fcn ( m, n, x, fvec, iflag )
            %      integer n
            %      real ( kind = rk ) fvec(m)
            %      integer iflag
            %      real ( kind = rk ) x(n)
            %
            %    N, is the number of variables. 
            %
            %    X(N), the point where the jacobian is evaluated. (p0)
            %
            %    FVEC(M), the functions evaluated at X.
            %
            %    LDFJAC, the leading dimension of FJAC, which must not be 
            %    less than M.
            %
            %  Output:
            %
            %    FJAC(LDFJAC,N), the M by N approximate jacobian matrix.
            %            
            epsfcn = 2.220446049250313e-16;
            epsmch = eps;
            epsilon = sqrt(max(epsfcn,epsmch));
            fjac = zeros([length(fvec),length(x)]);

            for j = 1:length(x)
                temp = x(j);
                h = epsilon * abs(temp);
                if h == 0.0 
                    h = epsilon;
                end

                x(j) = temp + h;
                args = [1, x(1), 0, x(2)];
                wa = fcn(x_data, y_data, args);
                
                x(j) = temp;
                fjac(:,j) = ( wa - fvec ) ./ h;
            end

        end

        function [aout, ipvt, rdiag, acnorm] = qrfac(a,pivot)
            %  QRFAC computes a QR factorization using Householder 
            %  transformations.
            %
            %  Discussion:
            %
            %    This function uses Householder transformations with 
            %    optional column pivoting to compute a QR factorization of 
            %    the M by N matrix A.  That is, QRFAC determines an 
            %    orthogonal matrix Q, a permutation matrix P, and an upper 
            %    trapezoidal matrix R with diagonal elements of 
            %    nonincreasing magnitude, such that A*P = Q*R.  
            %
            %    The Householder transformation for column 
            %    K, K = 1,2,...,min(M,N), is of the form
            %
            %      I - ( 1 / U(K) ) * U * U'
            %
            %    where U has zeros in the first K-1 positions.
            %
            %  Parameters:
            %
            %    M, the number of rows of A.
            %
            %    N, the number of columns of A.
            %
            %    A(LDA,N), the M by N array.
            %    A contains the matrix for which the QR factorization is to
            %    be computed.  
            %
            %    LDA, the leading dimension of A, which must
            %    be no less than M.
            %
            %    PIVOT, is TRUE if column pivoting is to be carried out.
            %
            %    LIPVT, the dimension of IPVT, which should 
            %    be N if pivoting is used.
            %
            %  Output:
            %
            %    aout, the strict upper trapezoidal part of A contains
            %    the strict upper trapezoidal part of R, and the lower 
            %    trapezoidal part of A contains a factored form of Q, 
            %    the non-trivial elements of the U vectors described above.
            %
            %    IPVT(LIPVT), defines the permutation matrix P 
            %    such that A*P = Q*R. Column J of P is column IPVT(J) of 
            %    the identity matrix. If PIVOT is false, IPVT is not 
            %    referenced.
            %
            %    RDIAG(N), contains the diagonal elements of R.
            %
            %    ACNORM(N), the norms of the corresponding
            %    columns of the input matrix A.  If this information is not 
            %    needed, then ACNORM can coincide with RDIAG.

            epsmch = eps;
            [m,n] = size(a);
            ipvt = zeros([1,n]);
            rdiag = zeros([1,n]);
            acnorm = zeros([1,n]);
            wa = zeros([1,n]);

            % Compute the initial column norms and initialize several 
            % arrays.
            for j = 1:n
                acnorm(j) = norm(a(1:m,j));
                rdiag(j) = acnorm(j);
                wa(j) = rdiag(j);
                if pivot               
                  ipvt(j) = j;
                end    
            end
            % Reduce A to R with Householder transformations.
            minmn = min( m, n );
            for j = 1:minmn
                %  Bring the column of largest norm into the pivot position
                if pivot
                    kmax = j;
                    for k = j:n
                        if rdiag(kmax) < rdiag(k)
                            kmax = k;
                        end
                    end

                    if kmax ~= j
                        r8_temp(1:m) = a(1:m,j);
                        a(1:m,j)     = a(1:m,kmax);
                        a(1:m,kmax)  = r8_temp(1:m);
                
                        rdiag(kmax) = rdiag(j);
                        wa(kmax) = wa(j);
                
                        i4_temp    = ipvt(j);
                        ipvt(j)    = ipvt(kmax);
                        ipvt(kmax) = i4_temp;
                    end
                end

                %  Compute the Householder transformation to reduce the
                %  J-th column of A to a multiple of the J-th unit vector.
                ajnorm = norm(a(:,j)); 
                if ajnorm ~= 0
                    if a(j,j) < 0
                        ajnorm = -ajnorm;
                    end
                    
                    a(j:m,j) = a(j:m,j)/ajnorm;
                    a(j,j) = a(j,j) + 1.0;

                    % Apply the transformation to the remaining columns and 
                    % update the norms.
                    for k = j+1:n
                        temp = dot(a(j:m,j), a(j:m,k) ) / a(j,j);
                        a(j:m,k) = a(j:m,k) - temp * a(j:m,j);
                        if pivot && rdiag(k) ~= 0
                            temp = a(j,k) / rdiag(k);
                            rdiag(k) = rdiag(k) * sqrt(max(0,1.0-temp^2));
                            if 0.05*(rdiag(k)/wa(k))^2 <= epsmch
                                rdiag(k) = sqrt(sum(a(j+1:m-j,k).^2));  
                                wa(k) = rdiag(k);
                            end
                        end
                    end
                end
                rdiag(j) = -ajnorm;
            end
            aout = a;
        end

        function [r, par, x, sdiag] = lmpar(r, ipvt, diag, qtb, delta, par)
            % LMPAR computes a parameter for the Levenberg-Marquardt method
            %
            %  Discussion:
            %
            %    Given an M by N matrix A, an N by N nonsingular diagonal
            %    matrix D, an M-vector B, and a positive number DELTA,
            %    the problem is to determine a value for the parameter
            %    PAR such that if X solves the system
            %
            %      A*X = B,
            %      sqrt ( PAR ) * D * X = 0,
            %
            %    in the least squares sense, and DXNORM is the euclidean
            %    norm of D*X, then either PAR is zero and
            %
            %      ( DXNORM - DELTA ) <= 0.1 * DELTA,
            %
            %    or PAR is positive and
            %
            %      abs ( DXNORM - DELTA) <= 0.1 * DELTA.
            %
            %    This function completes the solution of the problem
            %    if it is provided with the necessary information from the
            %    QR factorization, with column pivoting, of A.  That is, if
            %    A*P = Q*R, where P is a permutation matrix, Q has 
            %    orthogonal columns, and R is an upper triangular matrix 
            %    with diagonal elements of nonincreasing magnitude, then 
            %    LMPAR expects the full upper triangle of R, the 
            %    permutation matrix P, and the first N components of Q'*B.  
            %    On output LMPAR also provides an upper triangular matrix S 
            %    such that
            %
            %      P' * ( A' * A + PAR * D * D ) * P = S'* S.
            %
            %    S is employed within LMPAR and may be of separate interest
            %
            %    Only a few iterations are generally needed for convergence
            %    of the algorithm.  
            %
            %    If, however, the limit of 10 iterations is reached, then 
            %    the output PAR will contain the best value obtained so far
            %
            %  Parameters:
            %
            %    N, the order of R.
            %
            %    R(LDR,N),the N by N matrix.  The full upper triangle must 
            %    contain the full upper triangle of the matrix R.
            %
            %    LDR, the leading dimension of R.  LDR must be no less than
            %    N.
            %
            %    IPVT(N), defines the permutation matrix P such that 
            %    A*P = Q*R.  Column J of P is column IPVT(J) of the 
            %    identity matrix.
            %
            %    DIAG(N), the diagonal elements of the matrix D.
            %
            %    QTB(N), the first N elements of the vector Q'*B.
            %
            %    DELTA, an upper bound on the euclidean norm of D*X.  
            %      - DELTA should be positive.
            %
            %    PAR.  An initial estimate of the Levenberg-Marquardt 
            %    parameter.  PAR should be nonnegative.            
            %
            %  Output:
            %
            %    R(LDR,N) the full upper triangle is unaltered, and the 
            %    strict lower triangle contains the strict upper triangle 
            %    (transposed) of the upper triangular matrix S.
            %
            %    PAR.  The final estimate of the Levenberg-Marquardt 
            %    parameter. 
            %
            %    X(N), the least squares solution of the system
            %    A*X = B, sqrt(PAR)*D*X = 0, for the output value of PAR.
            %
            %    SDIAG(N), the diagonal elements of the upper triangular 
            %    matrix S.
    
            dwarf = realmin('single');
            [ldr,n] = size(r);
            x = zeros([1,n]);
            sdiag = zeros([1,n]);
            wa1 = zeros([1,n]);
            wa2 = zeros([1,n]);

            % Compute and store in X the Gauss-Newton direction.

            % If the jacobian is rank-deficient, obtain a least squares 
            % solution.
            nsing = n;

            for j = 1:n
                wa1(j) = qtb(j);
                if r(j,j) == 0.0 && nsing == n
                    nsing = j - 1;
                end
                if nsing < n
                    wa1(j) = 0.0;
                end
            end

            for k = 1:nsing
                j = nsing - k + 1;
                wa1(j) = wa1(j) / r(j,j);
                temp = wa1(j);
                wa1(1:j-1) = wa1(1:j-1) - r(1:j-1,j)*temp;
            end

            for j = 1:n
                l = ipvt(j);
                x(l) = wa1(j);
            end

            % Initialize the iteration counter.
            % Evaluate the function at the origin, and test
            % for acceptance of the Gauss-Newton direction.
            iter = 0;
            wa2(1:n) = diag(1:n) .* x(1:n);
            dxnorm = norm(wa2);
            fp = dxnorm - delta;

            if fp <= 0.1*delta
                if iter == 0
                    par = 0.0;
                end
                return
            end

            % If the jacobian is not rank deficient, the Newton
            % step provides a lower bound, PARL, for the zero of
            % the function.

            % Otherwise set this bound to zero.

            parl = 0.0;
            if n <= nsing 
                for j = 1:n
                    l = ipvt(j);
                    wa1(j) = diag(l) * ( wa2(l) / dxnorm );
                end

                for j = 1:n
                    sum2 = dot(wa1(1:j-1), r(1:j-1,j));
                    wa1(j) = ( wa1(j) - sum2 ) / r(j,j);
                end

                temp = norm(wa1);
                parl = ((fp/delta)/temp)/temp;
            end

            % Calculate an upper bound, PARU, for the zero of the function.
            for j = 1:n
                sum2 = dot( qtb(1:j), r(1:j,j) );
                l = ipvt(j);
                wa1(j) = sum2 / diag(l);
            end

            gnorm = enorm ( n, wa1 );
            paru = gnorm / delta;

            if paru == 0.0
                paru = dwarf/min(delta,0.1);
            end

            % If the input PAR lies outside of the interval (PARL, PARU),
            % set PAR to the closer endpoint.
            par = max(par,parl);
            par = min(par, paru);
            if par == 0.0
                par = gnorm / dxnorm;
            end

            % Beginning of iteration
            while true
                iter = iter + 1;
                % Evaluate the function at the current value of PAR.
                if par == 0.0
                    par = max(dwarf, 0.001*paru);
                end
                
                wa1(1:n) = sqrt(par) * diag(1:n);
                [r, x, sdiag] = ...
                    peaks_distribution.qrsolv(r, ipvt, diag, qtb);

                wa2 = diag * x;
                dxnorm = norm(wa2);
                temp = fp;
                fp = dxnorm - delta;

                % If the function is small enough, accept the current value 
                % of PAR.


            end
        end

        function [r, x, sdiag] = qrsolv(r, ipvt, diag, qtb)
            % QRSOLV solves a rectangular linear system A*x=b in the least 
            % squares sense.
            %
            %  Discussion:
            %
            %    Given an M by N matrix A, an N by N diagonal matrix D,
            %    and an M-vector B, the problem is to determine an X which
            %    solves the system
            %
            %      A*X = B
            %      D*X = 0
            %
            %    in the least squares sense.
            %
            %    This function completes the solution of the problem
            %    if it is provided with the necessary information from the
            %    QR factorization, with column pivoting, of A.  That is, if
            %    A*P = Q*R, where P is a permutation matrix, Q has 
            %    orthogonal columns, and R is an upper triangular matrix 
            %    with diagonal elements of nonincreasing magnitude, then 
            %    QRSOLV expects the full upper triangle of R, the 
            %    permutation matrix p, and the first N components of Q'*B.
            %
            %    The system is then equivalent to
            %
            %      R*Z = Q'*B
            %      P'*D*P*Z = 0
            %
            %    where X = P*Z.  If this system does not have full rank,
            %    then a least squares solution is obtained. On output 
            %    QRSOLV also provides an upper triangular matrix S such 
            %    that
            %
            %      P'*(A'*A + D*D)*P = S'*S.
            %
            %    S is computed within QRSOLV and may be of separate 
            %    interest.
            %
            %  Parameters:
            %
            %    R(LDR,N), the N by N matrix. The full upper triangle must 
            %    contain the full upper triangle of the matrix R.          
            %
            %    IPVT(N), defines the permutation matrix P such that 
            %      - A*P = Q*R.  
            %    Column J of P is column IPVT(J) of the identity matrix.
            %
            %    DIAG(N), the diagonal elements of the matrix D.
            %
            %    QTB(N), the first N elements of the vector Q'*B.            
            %
            %  Output:
            %
            %    R(LDR,N), the N by N matrix. The full upper triangle is 
            %    unaltered, and the strict lower triangle contains the 
            %    strict upper triangle (transposed) of the upper triangular 
            %    matrix S.
            %
            %    X(N), the least squares solution.
            %
            %    SDIAG(N), the diagonal elements of the upper triangular 
            %    matrix S.
            %

            [ldr,n] = size(r);
            x = zeros([1,n]);
            sdiag = zeros([1,n]);

            % Copy R and Q'*B to preserve input and initialize S.
            %
            % In particular, save the diagonal elements of R in X.
            for j = 1:n
                r(j:n) = r(j,j:n);
                x(j) = r(j,j);
            end

            wa = qtb;

            % Eliminate the diagonal matrix D using a Givens rotation.
            for j = 1:n
                % Prepare the row of D to be eliminated, locating the
                % diagonal element using P from the QR factorization.
                l = ipvt(j);
                if diag(l) ~= 0.0
                    sdiag(j:n) = 0.0;
                    sdiag(j) = diag(l);
                    % The transformations to eliminate the row of D
                    % modify only a single element of Q'*B
                    % beyond the first N, which is initially zero.
                    qtbpj = 0.0;
                    for k = j:n
                        % Determine a Givens rotation which eliminates the
                        % appropriate element in the current row of D.
                        if sdiag(k) ~= 0.0
                            if abs(r(k,k)) < abs(sdaig(k))
                                cotan = r(k,k)/sdiag(k);
                                s = 0.5/sqrt(0.25 + 0.25*cotan^2);
                                c = s * cotan;
                            else
                                t = sdaig(k) / r(k,k);
                                c = 0.5 / sqrt(0.25 + 0.25*t^2);
                                s = c * t;
                            end

                            % Compute the modified diagonal element of R 
                            % and the modified element of (Q'*B,0).
                            r(k,k) = c * r(k,k) + s * sdiag(k);
                            temp = c * wa(k) + s * qtbpj;
                            qtbpj = -s * wa(k) + c * qtbpj;
                            wa(k) = temp;

                            % Accumulate the tranformation in the row of S.
                            for i = k + 1:n
                                temp = c * r(i,k) + s * sdiag(i);
                                sdiag(i) = - s * r(i,k) + c * sdiag(i);
                                r(i,k) = temp;
                            end
                        end
                    end
                end

                % Store the diagonal element of S and restore the 
                % corresponding diagonal element of R.
                sdiag(j) = r(j,j);
                r(j,j) = x(j);
            end

            % Solve the triangular system for Z. If the system is 
            % singular, then obtain a least squares solution.
            nsing = n;
            for j = 1:n
                if sdiag(j) == 0.0 && nsing == n
                    nsing = j - 1;
                end

                if nsing < n
                    wa(j) = 0.0;
                end
            end

            for j = nsing:-1:1
                sum2 = dot(wa(j+1:nsing), r(j+1:nsing,j));
                wa(j) = (wa(j) - sum2)/ sdiag(j);
            end

            % Permute the components of Z back to components of X
            for j = 1:n
                l = ipvt(j);
                x(l) = wa(j);
            end
        end 
        
        function out = support_mask(x, arg)
            % This is just for the peaks over threshold method
            if arg < 0
                b = -1/arg;
            else
                b = inf;
            end
            a = 0;
            out = (a < x) & (x < b);
        end

    end
end


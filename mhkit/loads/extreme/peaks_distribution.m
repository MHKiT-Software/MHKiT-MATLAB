classdef peaks_distribution < handle  
    %PEAKS_DISTRIBUTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        qoi_peaks
        method
        paramEsts
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
            % Probability distribution of the peaks.

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
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
        end

        function peaks_over_threshold(obj, qoi_peaks)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
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
                ME = MException('MATLAB:extreme:ste_block_maxima:ppf',...
                    ['You must initialize the ste_block_maxima class with ' ...
                    'one of the available fit methods.\nMethods include\n' ...
                    '\tgev\n\tgumbel\n']);
                throw(ME);
            end
            
            [args, loc, scale] = obj.als();

            cond0 = scale > 0;
            cond1 = (0 < q) & (q < 1);
            cond = cond0 & cond1;

            if any(cond)
                a = args(1); c = args(2);
                out = -1*log(-q.^(1.0/a)+1).^(1.0/c);
                out = out * scale + loc;
            end

            if isempty(out)
                out = [];
            end
        end


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

            if isinf(loc_hat)
                loc_hat = 0;
            end
            if ~(~isinf(scale_hat) && 0 < scale_hat)
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

            [x0, ~, ~] = obj.als();            
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
                theta = ...
                    peaks_distribution.restore([1, 1.0, 0, 1], sim(k,:));
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
                    peaks_distribution.restore([1, 1.0, 0, 1], xr);
                fxr = obj.penalize_nnlf(theta);
                doshrink = false;

                if fxr < fsim(1)
                    xe = (1 + rho * chi) * xbar - rho * chi * sim(end,:);
                    theta = ...
                    peaks_distribution.restore([1, 1.0, 0, 1], xe);
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
                            theta = ...
                            peaks_distribution.restore([1, 1.0, 0, 1], xc);
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
                            peaks_distribution.restore([1, 1, 0, 1], xcc);
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
                                    [1, 1, 0, 1], sim(j,:));
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
            negxc = -x.^args(2);
            exm1c = -1*(exp(negxc)-1);
            logpdf = (log(args(1)) + log(args(2)) + ...
                (args(1)-1.0)*log(exm1c)...
                + negxc + (args(2)-1.0)*log(x));
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

        %% Utility Functions
        function [mu, mu2] = stats(obj)
            % statistics of the given RV using mean variance.            
            [pos, ~,~] = obj.als;
            loc = 0;
            scale = 1; 
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

        function [args, loc, scale] = als(obj)
            if strcmp(obj.method,'weibull')
                args = [obj.paramEsts{1},obj.paramEsts{2}];
                loc = obj.paramEsts{3};
                scale = obj.paramEsts{4};
            end
        end

    end

    methods(Static)
        function new_theta = restore(args, theta)
            i = 1;
            fixedn = [1, 3];
            for qq = 1:length(args)
                if ~ismember(qq,fixedn)
                    args(qq) = theta(i);
                    i = i + 1;
                end
            end
            new_theta = args;
        end
        
    end
end


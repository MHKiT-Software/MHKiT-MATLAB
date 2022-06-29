classdef genextreme < handle
    %GENEXTREME A generalized extreme value 
    % continuous random variable.
    %
    % Original function from scipy.stats._distn_infrastructure
    % Author:  Travis Oliphant  2002-2011 with contributions from
    %          SciPy Developers 2004-2011

    
    properties
        block_maxima
        paramEsts
        fcalls
        fit_called
    end
    
    methods
        function obj = genextreme(block_maxima)
            %GENEXTREME Construct an instance of this class
            %   Block maxima is an array of the maximum of each block being
            %   used to estimate the short-term extreme distribution
            obj.block_maxima = block_maxima;
            obj.paramEsts = cell([1,3]);
            obj.fcalls = 0;
            obj.fit_called = false;
        end

        function out = ppf(obj,q)
            % Percent point function (inverse of `cdf`) at q of the given 
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

            if ~obj.fit_called
                obj.fit();
            end
            
            args  = obj.paramEsts{1};
            loc   = obj.paramEsts{2};
            scale = obj.paramEsts{3};            
            cond0 = scale > 0;
            cond1 = (0 < q) && (q < 1);
            cond = cond0 && cond1;

            if cond
                x = -log(-log(q));
                if args ~= 0
                    out = -1*(exp(-args * x) - 1) / args;
                else
                    out = x;
                end
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

            args  = obj.paramEsts{1};
            loc   = obj.paramEsts{2};
            scale = obj.paramEsts{3};
            x = (x-loc)/scale;
            cond0 = scale > 0;
            cond1 = genextreme.support_mask(x, args) & scale > 0;
            cond = cond0 & cond1;
            out = zeros(size(cond));

            if any(cond)
                goodargs = x(cond);
                cx = goodargs*args;
                logex2 = log(1 + -cx);
                if args ~= 0
                    logpex2 = log(1 + -args*goodargs)/args;
                else
                    logpex2 = -goodargs;
                end
                pex2 = exp(logpex2);

                logpdf = zeros(size(goodargs));
                logpdf(cx==1 | cx == -inf) = -inf;
                logpdf(cx~=1 & cx ~= -inf) = -pex2 + logpex2 - logex2; 
                out(cond) = exp(logpdf)/scale;
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

            if ~obj.fit_called
                obj.fit();
            end

            args  = obj.paramEsts{1};
            loc   = obj.paramEsts{2};
            scale = obj.paramEsts{3};
            
            [a,b] = genextreme.get_support(args);
            lb = loc + a * scale;
            ub = loc + b * scale;

            out = obj.quad(lb,ub);

        end
        
        function fit(obj)
            % fit
            %   Return estimates of shape, location, and 
            %   scale parameters from data.
            %
            %   the fit is computed by minimizing the negative 
            %   log-likelihood function. A large, finite penalty
            %   (rather than infinite negative log-likelihood) is applied 
            %   for observations beyond the support of the distribution.

            if ~all(isfinite(obj.block_maxima))
                ME = MException('MATLAB:genextreme:fit',...
                    'The data contains non-finite values.');
                throw(ME);
            end

            % get distribution specific starting locations
            g = genextreme.skew(obj.block_maxima);
            if g < 0
                a = 0.5;
            else
                a = -0.5;
            end
            obj.paramEsts{1} = a;
            % Starting point for fit (shape arguments + loc + scale).
            [loc, scale] = obj.fit_loc_scale_support();
            obj.paramEsts{2} = loc;
            obj.paramEsts{3} = scale;
            [xopt, fopt, iter] = obj.fmin();
            obj.paramEsts{1} = xopt(1);
            obj.paramEsts{2} = xopt(2);
            obj.paramEsts{3} = xopt(3); 
            obj.fit_called = true;
        end   
        
        function [xopt, fopt, iter] = fmin(obj)
            % Minimize a function using the downhill simplex algorithm.
	        % 
            % This algorithm only uses function values, not derivatives or
            % second derivatives.
	        % 
            % Parameters
            % ----------        
	        % 
            % Returns
            % -------
            % xopt : ndarray
            %     Parameter that minimizes function.
            % fopt : float
            %     Value of function at minimum: ``fopt = func(xopt)``.
            % iter : int
            %     Number of iterations performed.    
	        %     
            % Notes
            % -----
            % Uses a Nelder-Mead simplex algorithm to find the minimum of 
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
	        %     
	        % 
            % References
            % ----------
            % .. [1] Nelder, J.A. and Mead, R. (1965), "A simplex method 
            %        for function minimization", The Computer Journal, 7,
            %        pp. 308-313
	        % 
            % .. [2] Wright, M.H. (1996), "Direct Search Methods: Once 
            %        Scorned, Now Respectable", in Numerical Analysis 1995, 
            %        Proceedings of the 1995 Dundee Biennial Conference in 
            %        Numerical Analysis, D.F. Griffiths and G.A. Watson 
            %        (Eds.), Addison Wesley Longman,Harlow, UK, pp. 191-208.
            xtol = 1e-4;
            ftol = 1e-4;
            obj.fcalls = 0;
            % initial guess
            x0 = [obj.paramEsts{1}, obj.paramEsts{2}, obj.paramEsts{3}];

            % Minimization of scalar function of one or more variables 
            % using the Nelder-Mead algorithm.
            rho = 1;
            chi = 2;
            psi = 0.5;
            sigma = 0.5;
            nonzdelt = 0.05;
            zdelt = 0.00025;
            N = length(x0);
            sim = zeros([N + 1, N]);
            sim(1,:) = x0;
            for k = 1:N
                y = x0;
                if y(k) ~= 0
                    y(k) = (1 + nonzdelt)*y(k);
                else
                    y(k) = zdelt;
                end
                sim(k + 1,:) = y;
            end
            maxiter = N * 200;
            maxfun = N * 200;
            one2np1 = 2:N + 1;
            fsim = zeros([N + 1,1]);

            for k= 1:(N + 1)
                fsim(k) = obj.penalized_nnlf(sim(k,:));                
            end
            [~,ind] = sort(fsim);
            fsim = fsim(ind);
            % sort so sim[0,:] has the lowest function value
            sim = sim(ind,:);
        
            iterations = 1;
            while obj.fcalls < maxfun && iterations < maxiter                
                if max(reshape(abs((sim(2:end,:) - sim(1,:))).',1,[]))...
                        <= 1e-4 && max(abs(fsim(1)-fsim(2:end))) <= 0.0001
                    break
                end

                xbar = sum(sim(1:end-1,:),1) / N;
                xr = (1 + rho) * xbar - rho * sim(end,:);
                fxr = obj.penalized_nnlf(xr);
                doshrink = false;

                if fxr < fsim(1)
                    xe = (1 + rho * chi) * xbar - rho * chi * sim(end,:);
                    fxe = obj.penalized_nnlf(xe);
                    if fxe < fxr
                        sim(end,:) = xe;
                        fsim(end,:) = fxe;
                    else
                        sim(end,:) = xr;
                        fsim(end,:) = fxr;
                    end
                else  % fsim(1) <= fxr
                    if fxr < fsim(end-1)
                        sim(end,:) = xr;
                        fsim(end,:) = fxr;
                    else % fxr >= fsim(end-1)
                    % Perform contraction
                        if fxr < fsim(end) 
                            xc = (1+psi * rho) * xbar-psi * rho*sim(end,:);
                            fxc = obj.penalized_nnlf(xc);
                        	if fxc <= fxr
                                sim(end,:) = xc;
                                fsim(end,:) = fxc;
                            else
                                doshrink = true;
                        	end
                        else
                            % Perform an inside contraction
                            xcc = (1 - psi) * xbar + psi * sim(end,:);
                            fxcc = obj.penalized_nnlf(xcc);
                            if fxcc < fsim(end)
                                sim(end,:) = xcc;
                                fsim(end,:) = fxcc;
                            else
                                doshrink = true;
                            end                                                   
                        end
                        if doshrink
                            for qq = 1:length(one2np1)
                                j = one2np1(j);
                                sim(j,:) = sim(1,:) + sigma * ...
                                    (sim(j,:) - sim(1,:));
                                fsim(j,:) = obj.penalized_nnlf(sim(j,:));
                            end
                        end
                    end
                end
                [~,ind] = sort(fsim);
                sim = sim(ind,:);
                fsim = fsim(ind,:);
                iterations = iterations + 1;
            end

            xopt = sim(1,:);
            fopt = min(fsim);
            iter = iterations;  
        end

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
            if obj.paramEsts{1} < 0
                a = 1.0 / min(obj.paramEsts{1}, -2.2250738585072014e-308);
                b = inf;
            else
                a = -inf;
                b = 1.0 / max(obj.paramEsts{1}, 2.2250738585072014e-308);
            end
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
            data_a = min(obj.block_maxima);
            data_b = max(obj.block_maxima);
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
                ME = MException('MATLAB:genextreme:fit_loc_scale_support',...
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

            muhat = mean(obj.block_maxima);
            mu2hat = var(obj.block_maxima);

            scale_hat = sqrt(mu2hat / mu2);
            loc_hat = muhat - scale_hat*mu ;

            if isinf(loc_hat)
                loc_hat = 0;
            end
            if ~(~isinf(scale_hat) && 0 < scale_hat)
                scale_hat = 1;
            end
        end

        function [mu, mu2] = stats(obj)
            % statistics of the given RV using mean variance.
            pos = obj.paramEsts{1};
            loc = 0;
            scale = 1;             
            
            [mu, mu2, ~, ~] = genextreme.stats2(pos);
            mu = mu * scale + loc;
            mu2 = mu2 * scale * scale;
        end
        
        function out = penalized_nnlf(obj, theta)
            % Penalized negative loglikelihood function.

            % i.e., - sum (log pdf(x, theta), axis=0) + penalty
            % where theta are the parameters (including loc and scale)
            obj.fcalls = obj.fcalls + 1;
            args = theta(1);
            loc = theta(2);
            scale = theta(3);
            x = ((obj.block_maxima - loc) ./ scale);
            n_log_scale = length(x) * log(scale);
            % return self._nnlf_and_penalty(x, args) + n_log_scale

            % Negative loglikelihood function
            cond0 = ~genextreme.support_mask(x, args);
            n_bad = sum(cond0(:));
            if n_bad > 0
                x = x(~cond0);
            end            
            cx = x*args;
            logex2 = log(1 + -1*cx);            
            if args ~= 0
                logpex2 = log(1 + -1*cx) /args;
            else
                logpex2 = -x;
            end
            pex2 = exp(logpex2);
            if args == 0
                logpex2(x == -inf) = 0.0;
            end
            logpdf = zeros(size(x));
            inds = (cx == 1) | (cx == -inf);
            sum_ = -pex2+logpex2-logex2;
            logpdf(inds) = -inf;
            logpdf(~inds) = sum_(~inds);
            if args == 1
                logpdf(x == 1) = 0.0;
            end
            finite_logpdf = isfinite(logpdf);
            n_bad = n_bad + sum(~finite_logpdf);
            if n_bad > 0
                penalty = n_bad * log(1.7976931348623157e+308) * 100.0;
                out = -sum(logpdf(finite_logpdf)) + penalty;
            else
                out = -sum(logpdf);
            end            
            % 
            out = out + n_log_scale;
        end       

        function out = quad(obj, a, b)
            % Compute a definite integral.
            % 
            % Integrate func from `a` to `b` (possibly infinite interval) 
            % using a technique from the Fortran library QUADPACK.
            
            args  = obj.paramEsts{1};

            % check the limits of integration: \int_a^b, expect a < b
            flip = b < a;
            a = min(a, b);
            b = max(a, b);
            infbounds = 0;

            if b ~= Inf && a ~= -Inf
                % standard integration
            elseif b == Inf && a ~= -Inf
                infbounds = 1;
                bound = a;
            elseif b == Inf && a == -Inf
                infbounds = 2;
                bound = 0;     % ignored
            elseif b ~= Inf && a == -Inf
                infbounds = -1;
                bound = b;
            else
                ME = MException('MATLAB:genextreme:quad',"Infinity " + ...
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

    methods(Static)

        function skewed = skew(data)
            % skew is third central moment / variance**(1.5)
            data = reshape(data.',1,[]);
            mu = mean(data);
            m2 = mean((data-mu).^2);
            m3 = mean((data-mu).^3);
            skewed = m3 / (m2^1.5);
        end

        function [m, v, sk, ku] = stats2(c)
            g1 = gamma(1*c + 1);
            g2 = gamma(2*c + 1);
            g3 = gamma(3*c + 1);
            g4 = gamma(4*c + 1);
            if abs(c) < 1e-7
                g2mg12 = (c*pi)^2./6.;
            else
                g2mg12 = g2-g1^2;
            end
            if abs(c) < 1e-7
                gam2k = pi^2./6.;
            else
                gam2k = (exp(log(abs(gamma(2.*c+1.))) - ...
                    2.*log(abs(gamma(c+1.)))) - 1)/c^2.;
            end
            eps = 1e-14;
            if abs(c) < eps
                gamk = 0.5772156649015329;
            else
                gamk = (exp( log(abs(gamma(c+1.))) ) - 1)/c;
            end

            if c < -1.0
                m = nan;
            else
                m = -gamk;
            end
            if c < -0.5
                v = nan;
            else
                v = g1^2 * gam2k;
            end

            % skewness
            if c >= -1.0/3.0
                sk1 = sign(c)*(-g3 + (g2 + 2*g2mg12)*g1)/g2mg12^1.5;
            else
                sk1 = nan;
            end
            if abs(c) <= eps^0.29
                sk = 12.*sqrt(6)*1.2020569031595942/pi^3;
            else
                sk = sk1;
            end
            
            % kurtosis
            if c >= -1.0/4.0
                ku1 = (g4 + (-4*g3 + 3*(g2 + g2mg12)*g1)*g1)/g2mg12^2;
            else
                ku1 = nan;
            end
            if abs(c) <= eps^0.23
                ku = 12.0/5.0;
            else
                ku = ku1-3.0;
            end
        end

        function out = support_mask(x, arg)
            if arg < 0
                a = 1.0/min(arg,-2.2250738585072014e-308);
                b = inf;
            else
                a = -inf;
                b = 1.0/max(arg, 2.2250738585072014e-308);
            end
            out = (a < x) & (x < b);
        end

        function [a,b] = get_support(arg)
            if arg < 0
                a = 1.0/min(arg,-2.2250738585072014e-308);
                b = inf;
            else
                a = -inf;
                b = 1.0/max(arg, 2.2250738585072014e-308);
            end
        end
    end
end


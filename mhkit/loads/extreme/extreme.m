classdef extreme < handle
    %EXTREME Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        block_maxima
        method
        paramEsts
        fcalls
    end
    
    methods
        function obj = extreme()
            obj.method = 'None';
        end

        %% Init methods for extreme class
        function ste_block_maxima_gev(obj, block_maxima)
            % Approximate the short-term extreme distribution using the 
            % block maxima method and the Generalized Extreme Value 
            % distribution.
            % 
            % Parameters
            % ----------
            % block_maxima: array
            %     Block maxima (i.e. largest peak in each block).
            % 
            % Returns
            % -------
            % Short-term extreme distribution paramters set in class
            % parameters
            %
            %   the fit is computed by minimizing the negative 
            %   log-likelihood function. A large, finite penalty
            %   (rather than infinite negative log-likelihood) is applied 
            %   for observations beyond the support of the distribution.

            obj.block_maxima = block_maxima;
            obj.method = 'general';
            obj.paramEsts = cell([1,3]);

            if ~all(isfinite(obj.block_maxima))
                ME = MException('MATLAB:extreme:ste_block_maxima_gev',...
                    'The data contains non-finite values.');
                throw(ME);
            end

            % get distribution specific starting locations
            g = extreme.skew(obj.block_maxima);
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
        end

        function ste_block_maxima_gumbel(obj, block_maxima)
            % Approximate the short-term extreme distribution using the 
            % block maxima method and the Gumbel (right) distribution.
            % 
            % Parameters
            % ----------
            % block_maxima: array
            %     Block maxima (i.e. largest peak in each block).
            % 
            % Returns
            % -------
            %   estimates of location and scale parameters from data.
            %
            %   The fit is computed by the method of maximum likelihood, 
            %   the estimators of the location and scale are the roots of 
            %   the equation defined in func and the value of the 
            %   expression for loc that follows. Source: Statistical 
            %   Distributions, 3rd Edition. Evans, Hastings, and Peacock 
            %   (2000), Page 101

            obj.block_maxima = block_maxima;
            obj.method = 'gumbel';
            obj.paramEsts = cell([1,2]);

            % Starting point for fit             
            [loc, scale] = obj.fit_loc_scale_support();
            obj.paramEsts{1} = loc;
            obj.paramEsts{2} = scale;
            [x, info] = obj.minpack_hybrd();  

            if info ~= 1
                if info == 0
                    ME = MException('MATLAB:extreme:minpack_hybrd',...
                    'Improper input parameters were entered.');
                    throw(ME);
                elseif info == 2
                    ME = MException('MATLAB:extreme:minpack_hybrd',...
                    ['The number of calls to function has reached ' ...
                    'maxfev = 400']);
                    throw(ME);
                elseif info == 3
                    ME = MException('MATLAB:extreme:minpack_hybrd',...
                    ['xtol=1e-12 is too small, no futher improvement in ' ...
                    'the approximate solution is possible']);
                    throw(ME);
                elseif info == 4
                    ME = MException('MATLAB:extreme:minpack_hybrd',...
                    ['The iteration is not making good progress, as ' ...
                    'measured by the \n  improvement from the last five ' ...
                    'Jacobian evaluations.']);
                    throw(ME);
                elseif info == 5
                    ME = MException('MATLAB:extreme:minpack_hybrd',...
                    ['The iteration is not making good progress, as ' ...
                    'measured by the \n  improvement from the last ten ' ...
                    'iterations.']);
                    throw(ME);
                end
            end
            
            obj.paramEsts{2} = x; % scale

            % Compute the log of the sum of exponentials of input elements
            a = -obj.block_maxima/x;
            a_max = max(a);
            a_max(isinf(a_max)) = 0;
            tmp = exp(a-a_max);
            s = sum(tmp);
            logsumexp = a_max + log(s);

            obj.paramEsts{1} = -x*(logsumexp - ...
                log(length(obj.block_maxima))); % loc
        end
        
        %% Statistical Methods
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

            if strcmp(obj.method, 'None')
                ME = MException('MATLAB:extreme:ppf',...
                    ['You must initialize the extreme class with one ' ...
                    'of the available fit methods.\nMethods include\n\' ...
                    'tste_block_maxima_gev\n\tste_block_maxima_gumbel\n']);
                throw(ME);
            end
            
            [args, loc, scale] = obj.als();

            cond0 = scale > 0;
            cond1 = (0 < q) && (q < 1);
            cond = cond0 && cond1;

            if any(cond)
                x = -log(-log(q));
                if isnan(args) || args == 0
                    out = x;
                else % args ~= 0
                    out = -1*(exp(-args * x) - 1) / args;
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

            if strcmp(obj.method, 'None')
                ME = MException('MATLAB:extreme:pdf',...
                    ['You must initialize the extreme class with one ' ...
                    'of the available fit methods.\nMethods include\n\' ...
                    'tste_block_maxima_gev\n\tste_block_maxima_gumbel\n']);
                throw(ME);
            end

            [args, loc, scale] = obj.als();

            x = (x-loc)/scale;
            cond0 = scale > 0;
            cond1 = extreme.support_mask(x, args) & scale > 0;
            cond = cond0 & cond1;
            out = zeros(size(cond));

            if any(cond)
                goodargs = x(cond);
                
                if isnan(args) || args == 0
                    logpex2 = -goodargs;
                    cx = zeros(size(goodargs));
                    logex2 = zeros(size(goodargs));
                else % args ~= 0
                    cx = goodargs*args;
                    logex2 = log(1 + -cx);
                    logpex2 = log(1 + -args*goodargs)/args;                
                end
                pex2 = exp(logpex2);

                logpdf = zeros(size(goodargs));
                logpdf(cx==1 | cx == -inf) = -inf;
                logpdf(cx~=1 & cx ~= -inf) = -pex2 + logpex2 - logex2; 
                out(cond) = exp(logpdf)/scale;
            end              
        end

        function out =  cdf(obj, x)
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
                ME = MException('MATLAB:extreme:cdf',...
                    ['You must initialize the extreme class with one ' ...
                    'of the available fit methods.\nMethods include\n\' ...
                    'tste_block_maxima_gev\n\tste_block_maxima_gumbel\n']);
                throw(ME);
            end

            [args, loc, scale] = obj.als();

            x = (x-loc)/scale;
            cond0 = scale > 0;
            cond1 = extreme.support_mask(x, args) & scale > 0;
            cond = cond0 & cond1;
            out = zeros(size(cond));

            if any(cond)
                goodargs = x(cond);
                if isnan(args) || args == 0
                    logpex2 = -goodargs;
                else % args ~= 0
                    logpex2 = log(1 + -args*goodargs)/args;                
                end

                out(cond) = exp(-exp(logpex2));
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
                ME = MException('MATLAB:extreme:expect',...
                    ['You must initialize the extreme class with one ' ...
                    'of the available fit methods.\nMethods include\n\' ...
                    'tste_block_maxima_gev\n\tste_block_maxima_gumbel\n']);
                throw(ME);
            end

            [args, loc, scale] = obj.als();
            
            [a,b] = extreme.get_support(args);
            lb = loc + a * scale;
            ub = loc + b * scale;

            out = obj.quad(lb,ub);
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
            [args, ~,~] = obj.als;
            [a, b] = extreme.get_support(args);
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
                ME = MException('MATLAB:extreme:fit_loc_scale_support',...
                    'Runtime Error');
                throw(ME);
            end

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
    
        function [retval, info] = minpack_hybrd(obj)
            % By the method of maximum likelihood, the estimators of the
            % location and scale are the roots of the equation defined in
            % func and the value of the expression for loc that follows.
            % Source: Statistical Distributions, 3rd Edition. Evans, 
            % Hastings, and Peacock (2000), Page 101
            %  
            % Find the roots of a multivariate function using hybrd 
            % routine (modified Powell method).
            % https://people.sc.fsu.edu/~jburkardt/f_src/minpack/minpack.f90
            %
            %  HYBRD seeks a zero of N nonlinear equations in N variables.
            %
            %  Discussion:
            %
            %    HYBRD finds a zero of a system of N nonlinear functions in 
            %    N variables by a modification of the Powell hybrid method.  
            %    The user must provide a subroutine which calculates the 
            %    functions.            
            %    The jacobian is then calculated by a forward-difference 
            %    approximation.
            %
            % func: the name of the user-supplied subroutine which 
            %      calculates the functions.
            % n: the number of functions and variables
            % x: an initial estimate of the solution vector.
            % XTOL: Termination occurs when the relative error between 
            %       two consecutive iterates is at most XTOL
            % MAXFEV: Termination occurs when the number of calls to FCN 
            %         is at least MAXFEV
            % ML, MU: specify the number of subdiagonals and superdiagonals 
            %         within the band of the jacobian matrix
            % EPSFCN: is used in determining a suitable step length for the 
            %         forward-difference approximation.
            
            %func = @extreme.max_like; 
            retval = nan;
            x = obj.paramEsts{2}; % x0 initial value
            n = length(x);
            xtol = 1e-12;
            maxfev = 200 * (n + 1);
            ml = n-1;
            mu = n-1;
            mode = 1;
            fact = 100;
            diag = zeros([n,1]);
            wa3 = zeros([n,1]);
            r = zeros([(n*(n+1))/2,1]);
            epsfcn = 2.220446049250313e-16;
            
            %  Evaluate the function at the starting point
            %  and calculate its norm.
            fvec = extreme.max_like(x,obj.block_maxima);
            nfev = 1;
            fnorm = sqrt(sum(n^2));

            % Determine the number of calls to FCN needed to compute the 
            % jacobian matrix.
            msum = min( ml + mu + 1, n );
            
            %  Initialize iteration counter and monitors.
            info = 0;
            iter = 1;
            ncsuc = 0;
            ncfail = 0;
            nslow1 = 0;
            nslow2 = 0;
            
            % Beginning of the outer loop.
            while true
                jeval = true;
                % Calculate the jacobian matrix.
                fjac = obj.fdjac1(n, x, fvec, n, ml, mu);
                nfev = nfev + msum;
    
                % Compute the QR factorization of the jacobian.                
                [fjac, iwa, wa1, wa2] = ...
                    extreme.qrfac( n, n, fjac, false, 1);

                %  On the first iteration, if MODE is 1, scale according
                %  to the norms of the columns of the initial jacobian.
                if iter == 1
                    if mode ~= 2
                        diag(1:n) = wa2(1:n);
                        for j = 1:n
                            if wa2(j) == 0
                                diag(j) = 1.0;
                            end
                        end
                    end

                    %  On the first iteration, calculate the norm of the 
                    %  scaled X and initialize the step bound DELTA.
                    wa3(1:n) = diag(1:n) * x(1:n);
                    xnorm = norm(wa3);
                    delta = fact * xnorm;
                    if delta == 0.0
                        delta = fact;
                    end 
                end % if iter=1

                % Form Q' * FVEC and store in QTF.
                qtf(1:n) = fvec(1:n);
                for j = 1:n
                    if fjac(j,j) ~= 0.0
                        temp = -1*dot(qtf(j:n), fjac(j:n,j) ) / fjac(j,j);
                        qtf(j:n) = qtf(j:n) + fjac(j:n,j) * temp;
                    end
                end

                % Copy the triangular factor of the QR factorization into R
                sing = false;

                for j = 1:n
                    l = j;
                    for i = 1:j-1
                        r(l) = fjac(i,j);
                        l = l + n - i;
                    end
                    r(l) = wa1(j);
                    if wa1(j) == 0.0
                        sing = true;
                    end
                end

                % Accumulate the orthogonal factor in FJAC.
                fjac = extreme.qform(n, n, fjac);

                % Rescale if necessary.
                if mode ~= 2
                    for j = 1:n
                        diag(j) = max(diag(j),wa2(j));
                    end
                end

                % Beginning of the inner loop.
                while true
                    % Determine the direction P.
                    wa1 = extreme.dogleg(n, r, diag, qtf, delta);

                    % Store the direction P and X + P.
                    % Calculate the norm of P.
                    wa1(1:n) = -wa1(1:n);
                    wa2(1:n) = x(1:n) + wa1(1:n);
                    wa3(1:n) = diag(1:n) * wa1(1:n);
                    pnorm = norm(wa3);

                    % On the first iteration, adjust the initial step bound
                    if iter == 1
                        delta = min(delta, pnorm);
                    end

                    % Evaluate the function at X + P and calculate its norm
                    wa4 = extreme.max_like(wa2,obj.block_maxima);
                    nfev = nfev + 1;
                    fnorm1 = norm(wa4);

                    % Compute the scaled actual reduction.
                    actred = -1.0;
                    if fnorm1 < fnorm
                        actred = 1.0 - (fnorm1/fnorm)^2;
                    end

                    % Compute the scaled predicted reduction.
                    l = 1;
                    for i = 1:n
                        sum2 = 0.0;
                        for j = 1:n
                            sum2 = sum2 + r(l)*wa1(j);
                            l = l + 1;
                        end
                        wa3(i) = qtf(i) + sum2;
                    end

                    temp = norm(wa3);
                    prered = 0.0;
                    if temp < fnorm 
                        prered = 1.0 - (temp/fnorm)^2;
                    end

                    % Compute the ratio of the actual to the predicted 
                    % reduction.
                    if prered > 0.0
                        ratio = actred / prered;
                    else
                        ratio = 0;
                    end

                    % Update the step bound.
                    if ratio < 0.1
                        ncsuc = 0;
                        ncfail = ncfail + 1;
                        delta = 0.5*delta;
                    else
                        ncfail = 0;
                        ncsuc = ncsuc + 1;
                        if ratio >= 0.5 || ncsuc > 1
                            delta = max(delta,pnorm/0.5);
                        end
                        if abs(ratio-1.0) <= 0.1
                            delta = pnorm/0.5;
                        end
                    end

                    % On successful iteration, update X, FVEC, and their 
                    % norms.
                    if ratio >= 0.0001
                        x(1:n) = wa2(1:n);
                        wa2(1:n) = diag(1:n) * x(1:n);
                        fvec(1:n) = wa4(1:n);
                        xnorm = norm(wa2);
                        fnorm = fnorm1;
                        iter = iter + 1;
                    end

                    % Determine the progress of the iteration.
                    nslow1 = nslow1 + 1;
                    if actred >= 0.001
                        nslow1 = 0;
                    end

                    if jeval
                        nslow2 = nslow2 + 1;
                    end

                    if actred >= 0.1
                        nslow2 = 0;
                    end

                    % Test for convergence 
                    if delta <= xtol*xnorm || fnorm == 0.0
                        info = 1;
                        retval = x;
                        return
                    end

                    % Tests for termination and stringent tolerances.
                    if maxfev <= nfev
                        info = 2;
                        return
                    end
                    if (0.1*max(0.1*delta,pnorm)) <= eps*xnorm
                        info = 3;
                        return
                    end
                    if nslow2 == 5
                        info = 4;
                        return
                    end
                    if nslow1 == 10
                        info = 5;
                        return
                    end

                    % Criterion for recalculating jacobian.
                    if ncfail == 2
                        break
                    end

                    % Calculate the rank one modification to the jacobian
                    % and update QTF if necessary.
                    for j = 1:n
                        sum2 = dot(wa4(1:n),fjac(1:n,j));
                        wa2(j) = (sum2 - wa3(j))/pnorm;
                        wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm);
                        if ratio >= 0.0001
                            qtf(j) = sum2;
                        end
                    end

                    % Compute the QR factorization of the updated jacobian                    
                    [r, wa2, wa3, sing] =...
                        extreme.r1updt(n, n, r, wa1, wa2);                    
                    fjac = extreme.r1mpyq(n, n, fjac, wa2, wa3);
                    qtf = extreme.r1mpyq(1, n, qtf, wa2, wa3);
                
                    jeval = false; 

                end % while inner loop                
            end % while outer loop
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
            % Variance (matlab's var() function does not return the correct
	        % value
            mu2hat = mean(abs(obj.block_maxima -...
		        mean(obj.block_maxima)).^2);

            scale_hat = sqrt(mu2hat / mu2);
            loc_hat = muhat - scale_hat*mu ;

            if isinf(loc_hat)
                loc_hat = 0;
            end
            if ~(~isinf(scale_hat) && 0 < scale_hat)
                scale_hat = 1;
            end
        end
       
        function fjac = fdjac1(obj, n, x, fvec, ldfjac, ml, mu)
            %  FDJAC1 estimates an N by N jacobian matrix using forward 
            %  differences.
            %
            %  Discussion:
            %
            %    This function computes a forward-difference approximation
            %    to the N by N jacobian matrix associated with a specified
            %    problem of N functions in N variables. If the jacobian has
            %    a banded form, then function evaluations are saved by only
            %    approximating the nonzero terms.
            %
            %  Parameters:
            %
            %    N, the number of functions and variables.
            %
            %    X(N), the point where the jacobian is evaluated.
            %
            %    FVEC(N), the functions evaluated at X.                     
            %
            %    LDFJAC, the leading dimension of FJAC, which must not be 
            %    less than N.
            %
            %    ML, MU, specify the number of subdiagonals and
            %    superdiagonals within the band of the jacobian matrix.  
            %    If the jacobian is not banded, set ML and MU to N-1.
            %
            %  Output:
            %
            %    FJAC(LDFJAC,N), the N by N approximate jacobian matrix.
            
            epsfcn = 2.220446049250313e-16;
            epsmch = eps;
            epsilon = sqrt(max(epsfcn,epsmch));
            msum = ml + mu + 1;

            fjac = zeros(ldfjac,n);
            % Computation of dense approximate jacobian.
            if n <= msum
                for j = 1:n
                    temp = x(j);
                    h = epsilon * abs(temp);
                    if h == 0.0
                        h = epsilon;
                    end

                    x(j) = temp + h;
                    wa1 = extreme.max_like(x, obj.block_maxima);

                    x(j) = temp;
                    fjac(1:n,j) = (wa1(1:n) - fvec(1:n))./h;
                end            
            else % Computation of banded approximate jacobian
                for k = 1:msum
                    for j = k:msum:n
                        wa2(j) = x(j);
                        h = epsilon * abs(wa2(j));
                        if h == 0.0
                            h = epsilon;
                        end
                        x(j) = wa2(j) + h;
                    end

                    wa1 = extreme.max_like(x, obj.block_maxima);
                    for j = k:msum:n
                        x(j) = wa2(j);

                        h = epsilon * abs(wa2(j));
                        if h == 0.0
                            h = epsilon;
                        end

                        for i = 1:n
                            if (j-mu) <= i && i <= (j+ml)
                                fjac(i,j) = (wa1(i)-fvec(i))./h;
                            end
                        end
                    end
                end
            end
        end

        %% Utility Functions
        function [mu, mu2] = stats(obj)
            % statistics of the given RV using mean variance.
            sup_methods = {'general', 'gumbel'};
            if any(strcmp(sup_methods,obj.method))
                [pos, ~,~] = obj.als;
                loc = 0;
                scale = 1;             
                
                if isnan(pos)
                    mu = 0.5772156649015329;
	                mu2 = (pi*pi/6.0);
                else
                    [mu, mu2, ~, ~] = extreme.stats2(pos);
                end
                mu = mu * scale + loc;
                mu2 = mu2 * scale * scale;
            else
                ME = MException('MATLAB:extreme:stats',...
                    'stats not supported for method: %s', obj.method);
                throw(ME);
            end
        end

        function [args, loc, scale] = als(obj)
            if strcmp(obj.method,'general')
                args = obj.paramEsts{1};
                loc = obj.paramEsts{2};
                scale = obj.paramEsts{3};
            elseif strcmp(obj.method,'gumbel')
                args = nan;
                loc = obj.paramEsts{1};
                scale = obj.paramEsts{2};
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
                ME = MException('MATLAB:extreme:quad',"Infinity " + ...
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
            if isnan(arg)
                a = -inf;
                b = inf;                
            elseif arg < 0
                a = 1.0/min(arg,-2.2250738585072014e-308);
                b = inf;
            else
                a = -inf;
                b = 1.0/max(arg, 2.2250738585072014e-308);
            end
            out = (a < x) & (x < b);
        end

        function [a,b] = get_support(arg)
            if isnan(arg)
                a = -inf;
	            b = inf;  
            elseif arg < 0
                a = 1.0/min(arg,-2.2250738585072014e-308);
                b = inf;
            else
                a = -inf;
                b = 1.0/max(arg, 2.2250738585072014e-308);
            end
        end
    
        function [aout, ipvt, rdiag, acnorm] = qrfac(m,n,a,pivot,lipvt)
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
            ipvt = zeros([lipvt,1]);
            rdiag = zeros([n,1]);
            acnorm = zeros([n,1]);
            wa = zeros([n,1]);

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
                ajnorm = norm (a(j,j));
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
                                rdiag(k) = norm( a(j+1,k));
                                wa(k) = rdiag(k);
                            end
                        end
                    end
                end
                rdiag(j) = -ajnorm;
            end
            aout = a;
        end

        function out = qform(m, n, q)
            %  Purpose:
            %
            %    qform() constructs the standard form of Q from its 
            %    factored form.
            %
            %  Discussion:
            %
            %    This function proceeds from the computed QR factorization of
            %    an M by N matrix A to accumulate the M by M orthogonal matrix
            %    Q from its factored form.            
            %
            %  Parameters:
            %
            %    m, the number of rows of A, and the order of Q.
            %
            %    n, the number of columns of A.
            %
            %    Q[LDQ*N], the full lower trapezoid in the first min(M,N) 
            %    columns of Q contains the factored form.            
            %
            %    LDQ, the leading dimension of the array Q.
            %
            %  Output:
            %
            %    out, Q that has been accumulated into a square matrix.

            % Zero out the upper triangle of Q in the first min(M,N) 
            % columns.
            minmn = min( m, n );
            wa = zeros([m,1]);
            out = zeros(size(q));

            for j = 2:minmn
                q(1:j-1,j) = 0.0;
            end

            % Initialize remaining columns to those of the identity matrix.
            q(1:m, n+1:m) = 0.0;

            for j = n+1:m
                q(j,j) = 1.0;
            end

            % Accumulate Q from its factored form.
            for l = 1:minmn
                k = minmn - l + 1;
                wa(k:m) = q(k:m,k);
                q(k:m,k) = 0.0;
                q(k,k) = 1.0;

                if wa(k) ~= 0.0
                    for j = k:m
                        temp = dot(wa(k:m), q(k:m,j) ) / wa(k);
                        q(k:m,j) = q(k:m,j) - temp * wa(k:m);
                    end
                end
            end

            out = q;
        end

        function x = dogleg(n, r, diag, qtb, delta)
            %  Purpose:
            %
            %    dogleg() combines Gauss-Newton and gradient for a 
            %    minimizing step.
            %
            %  Discussion:
            %
            %    Given an M by N matrix A, an n by n nonsingular diagonal
            %    matrix d, an m-vector b, and a positive number delta, the
            %    problem is to determine the convex combination x of the
            %    gauss-newton and scaled gradient directions that minimizes
            %    (a*x - b) in the least squares sense, subject to the
            %    restriction that the euclidean norm of d*x be at most 
            %    delta.
            %
            %    This function completes the solution of the problem
            %    if it is provided with the necessary information from the
            %    qr factorization of a. 
            %
            %    That is, if a = q*r, where q has orthogonal columns and r 
            %    is an upper triangular matrix, then dogleg expects the 
            %    full upper triangle of r and the first n components of 
            %    Q'*b.
            %
            %  Parameters:
            %
            %    n, the order of R.
            %
            %    r[lr], the upper triangular matrix R stored by rows.           
            %
            %    diag[n], the diagonal elements of the matrix D.
            %
            %    qtb[n], the first n elements of the vector 
            %    (q transpose)*b.
            %
            %    delta, an upper bound on the euclidean norm of d*x.
            %
            %  Output
            %
            %    x[n], contains the desired convex combination of the 
            %    gauss-newton direction and the scaled gradient direction.

            epsmch = eps;
            x = zeros([n,1]);
            wa1 = zeros([n,1]);
            wa2 = zeros([n,1]);

            % Calculate the Gauss-Newton direction.
            jj = ( n * ( n + 1 ) ) / 2 + 1;

            for k = 1:n
                j = n - k + 1;
                jj = jj - k;
                l = jj + 1;
                sum2 = 0.0;

                for i = j+1:n
                    sum2 = sum2 + r(l)*x(i);
                    l = l + 1;
                end

                temp = r(jj);

                if temp == 0.0
                    l = j;
                    for i = 1:j
                        temp = max(temp, abs(r(l)));
                        l = l + n - i;
                    end

                    if temp == 0.0
                        temp = epsmch;
                    else
                        temp = epsmch*temp;
                    end
                end

                x(j) = ( qtb(j) - sum2 ) / temp;
            end

            % Test whether the Gauss-Newton direction is acceptable.
            wa1(1:n) = 0.0;
            wa2(1:n) = diag(1:n) * x(1:n);
            qnorm = norm(wa2);

            if qnorm <= delta
                return
            end

            % The Gauss-Newton direction is not acceptable.
            % Calculate the scaled gradient direction.
            l = 1;
            for j = 1:n
                temp = qtb(j);
                for i = j:n
                    wa1(i) = wa1(i) + r(l) * temp;
                    l = l + 1;
                end
                wa1(j) = wa1(j) / diag(j);
            end

            % Calculate the norm of the scaled gradient.
            % Test for the special case in which the scaled gradient is 
            % zero.
            gnorm = norm(wa1);
            sgnorm = 0.0;
            alpha = delta./qnorm;

            if gnorm ~= 0.0
                % Calculate the point along the scaled gradient which 
                % minimizes the quadratic.
                wa1(1:n) = ( wa1(1:n) / gnorm ) / diag(1:n);
                l = 1;
                for j = 1:n
                  sum2 = 0.0;
                  for i = j:n
                    sum2 = sum2 + r(l) * wa1(i);
                    l = l + 1;
                  end 
                  wa2(j) = sum2;
                end 
            
                temp = enorm ( n, wa2 );
                sgnorm = ( gnorm / temp ) / temp;
                alpha = 0.0;

                % If the scaled gradient direction is not acceptable,
                % calculate the point along the dogleg at which the
                % quadratic is minimized.
                if sgnorm < delta
                    bnorm = norm(qtb);
                    temp = (bnorm/gnorm) * (bnorm/qnorm) * (sgnorm/delta);
                    temp = temp - (delta/qnorm) * (sgnorm/delta)^2 ...
                        + sqrt((temp - (delta/qnorm) )^2 ...
                        + (1.0 - (delta/qnorm)^2)*(1.0 - (sgnorm/delta)^2));
                    alpha = ((delta/qnorm)*(1.0 - (sgnorm/delta)^2 ))/temp;
                end
            end

            % Form appropriate convex combination of the Gauss-Newton
            % direction and the scaled gradient direction.
            temp = (1.0 - alpha) * min(sgnorm, delta);
            x(1:n) = temp*wa1(1:n) + alpha*x(1:n);
        end

        function [sout, vout, w, sing] = r1updt(m, n, s, u, v)
            % R1UPDT re-triangularizes a matrix after a rank one update.
            %
            %  Discussion:
            %
            %    Given an M by N lower trapezoidal matrix S, an M-vector U, 
            %    and an N-vector V, the problem is to determine an 
            %    orthogonal matrix Q such that
            %
            %      (S + U * V' ) * Q
            %
            %    is again lower trapezoidal.
            %
            %    This function determines Q as the product of 2 * (N - 1)
            %    transformations
            %
            %      GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
            %
            %    where GV(I), GW(I) are Givens rotations in the (I,N) plane
            %    which eliminate elements in the I-th and N-th planes,
            %    respectively.  Q itself is not accumulated, rather the
            %    information to recover the GV and GW rotations is returned
            %
            %  Parameters:
            %
            %    m, the number of rows of S.
            %
            %    n, the number of columns of S. N must not exceed M.
            %
            %    s(ls), the lower trapezoidal matrix S stored by columns. 
            %
            %    u(m), the U vector.
            %
            %    v(n),  On input, V must contain the vector V.
            %
            %  Output:
            %
            %    sout contains the lower trapezoidal
            %    matrix produced as described above.
            %
            %    vout contains the information necessary to recover the
            %    Givens rotations GV described above.
            %
            %    w(m), contains information necessary to
            %    recover the Givens rotations GW described above.
            %
            %    SING, is set to TRUE if any of the diagonal elements
            %    of the output S are zero.  Otherwise SING is set FALSE.

            % Initialize the diagonal element pointer.
            jj = (n*(2*m - n + 1)) / 2 - (m - n);

            % Move the nontrivial part of the last column of S into W.
            l = jj;
            for i = n:m
                w(i) = s(l);
                l = l + 1;
            end

            % Rotate the vector V into a multiple of the N-th unit vector
            % in such a way that a spike is introduced into W.
            for j = n-1:-1:1
                jj = jj - (m - j +1);
                w(j) = 0.0;
                if v(j) ~= 0.0
                    % Determine a Givens rotation which eliminates the
                    % J-th element of V.
                    if abs(v(n)) < abs(v(j))                        
                        cotan = v(n) / v(j);
                        sn = 0.5/sqrt(0.25 + 0.25*cotan^2);
                        cs = sn * cotan;
                        tau = 1.0;
                        if abs(cs)*inf > 1.0
                            tau = 1.0/cs;
                        end
                    else
                        tn = v(j)/v(n);
                        cs = 0.5 / sqrt(0.25+0.25*tn^2);
                        sn = cs*tn;
                        tau = sn;
                    end

                    % Apply the transformation to V and store the 
                    % information necessary to recover the Givens rotation
                    v(n) = sn * v(j) + cs * v(n);
                    v(j) = tau;

                    % Apply the transformation to S and extend the spike 
                    % in W.
                    l = jj;
                    for i = j:m
                        temp = cs * s(l) - sn*w(i);
                        w(i) = sn * s(l) + cs*w(i);
                        s(l) = temp;
                        l = l + 1;
                    end
                end
            end

            % Add the spike from the rank 1 update to W.
            w(1:m) = w(1:m) + v(n) * u(1:m);

            % Eliminate the spike.
            sing = false;

            for j = 1:n-1
                if w(j) ~= 0.0
                    % Determine a Givens rotation which eliminates the
                    % J-th element of the spike.
                    if abs(s(jj)) < abs(w(j))
                        cotan = s(jj) / w(j);
                        sn = 0.5 / sqrt(0.25+0.25*cotan^2);
                        cs = sn * cotan;
                        if 1.0 < abs(cs)*inf
                            tau = 1.0 / cs;
                        else
                            tau = 1.0;
                        end
                    else
                        tn = w(j) / s(jj);
                        cs = 0.5 / sqrt(0.25+0.25*tn^2);
                        sn = cs*tn;
                        tau = sn;
                    end

                    % Apply the transformation to S and reduce the spike 
                    % in W.
                    l = jj;
                    for i = j:m
                        temp = cs*s(l) + sn*w(i);
                        w(i) = -sn*s(l) + cs*w(i);
                        s(l) = temp;
                        l = l + 1;
                    end

                    % Store the information necessary to recover the 
                    % Givens rotation.
                    w(j) = tau;
                end

                % Test for zero diagonal elements in the output S
                if s(jj) == 0.0
                    sing = true;
                end

                jj = jj + (m-j+1);
            end

            % Move W back into the last column of the output S.
            l = jj;
            for i = n:m
                s(l) = w(i);
                l = l + 1;
            end

            if s(jj) == 0.0
                sing = true;
            end

            % assign output
            sout = s;
            vout = v;
        end

        function a_out = r1mpyq(m, n, a, v, w)
            % r1mpyq computes A*Q, where Q is the product of Householder 
            % transformations.
            %
            %  Discussion:
            %
            %    Given an M by N matrix A, this function computes A*Q where
            %    Q is the product of 2*(N - 1) transformations
            %
            %      GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
            %
            %    and GV(I), GW(I) are Givens rotations in the (I,N) plane which
            %    eliminate elements in the I-th and N-th planes, respectively.
            %    Q itself is not given, rather the information to recover the
            %    GV, GW rotations is supplied.
            %
            %  Parameters:
            %
            %    m, the number of rows of A.
            %
            %    n, the number of columns of A.
            %
            %    a(LDA,N), the M by N array. the matrix A to be 
            %    postmultiplied by the orthogonal matrix Q.            
            %            
            %    v(N), w(N), contain the information necessary
            %    to recover the Givens rotations GV and GW.
            %
            %  Output:
            %
            %    a_out, the value of A*Q.

            % Apply the first set of Givens rotations to A.
            for j = n-1:-1:1
                if abs(v(j)) > 1.0
                    c = 1.0/v(j);
                    s = sqrt(1.0-c^2);
                else
                    s = v(j);
                    c = sqrt(1.0-s^2);
                end

                for i = 1:m
                    temp = c*a(i,j) - s*a(i,n);
                    a(i,n) = s*a(i,j) + c*a(i,n);
                    a(i,j) = temp;
                end
            end

            % Apply the second set of Givens rotations to A.
            for j = 1:n-1
                if abs(w(j)) > 1.0
                    c = 1.0/w(j);
                    s = sqrt(1.0 - c^2);
                else
                    s = w(j);
                    c = sqrt(1.0 - s^2);
                end

                for i = 1:m
                    temp = c*a(i,j) + s*a(i,n);
                    a(i,n) = -s*a(i,j) + c*a(i,n);
                    a(i,j) = temp;
                end
            end

            a_out = a;
        end

        function out = max_like(scale, data)
            % method of maximum likelihood
            sdata = -data/scale;
            maxlogw = max(sdata);
            weights = exp(sdata - maxlogw);
            wavg = sum(data .* weights) / sum(weights);
            out = mean(data) - wavg - scale;
        end        
    end
end
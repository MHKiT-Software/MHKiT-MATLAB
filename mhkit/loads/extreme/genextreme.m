classdef genextreme
    %GENEXTREME A generalized extreme value 
    % continuous random variable.
    %
    % Original function from scipy.stats._distn_infrastructure
    % Author:  Travis Oliphant  2002-2011 with contributions from
    %          SciPy Developers 2004-2011

    
    properties
        block_maxima
        paramEsts
    end
    
    methods
        function obj = genextreme(block_maxima)
            %GENEXTREME Construct an instance of this class
            %   Block maxima is an array of the maximum of each block being
            %   used to estimate the short-term extreme distribution
            obj.block_maxima = block_maxima;
            paramEsts = cell([1,3]);
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

            start = nan([1,2]);
            % get distribution specific starting locations
            g = obj.skew(obj.block_maxima);
            if g < 0
                a = 0.5;
            else
                a = -0.5;
            end
            obj.paramEsts{1} = a;
            % Starting point for fit (shape arguments + loc + scale).
            [loc, scale] = obj.fit_loc_scale_support(obj.block_maxima);
            obj.paramEsts{2} = loc;
            obj.paramEsts{3} = scale;
        end   

        function skewed = skew(obj, data)
            % skew is third central moment / variance**(1.5)
            data = reshape(data.',1,[]);
            mu = mean(data);
            m2 = mean((data-mu).^2);
            m3 = mean((data-mu).^3);
            skewed = m3 / (m2^1.5);
        end

        function [loc_hat, scale_hat] = fit_loc_scale_support(obj, data)
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

            % Estimate location and scale according to the method of 
            % moments.
            [loc_hat, scale_hat] = obj.fit_loc_scale(data);

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
            data_a = min(data);
            data_b = max(data);
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

        function [loc_hat, scale_hat] = fit_loc_scale(obj, data)
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
            
            [mu, mu2, ~, ~] = obj.stats2(pos);
            mu = mu * scale + loc;
            mu2 = mu2 * scale * scale;
        end

        function [m, v, sk, ku] = stats2(obj, c)
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


    end
end


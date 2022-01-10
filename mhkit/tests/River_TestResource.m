classdef River_TestResource < matlab.unittest.TestCase

    methods (Test) 

        function test_Froude_number(testCase)
            v = 2;
            h = 5;
            Fr = Froude_number(v, h);
            assertEqual(testCase,Fr, 0.286,'AbsTol',0.001);
        end

        function test_exceedance_probability(testCase)
            % Create arbitrary discharge between 0 and 8(N=9)
            Q = struct('Discharge',[0;1;2;3;4;5;6;7;8],'time',[0 1 2 3 4 5 6 7 8]);
            % Rank order for non-repeating elements simply adds 1 to each element
            %if N=9, max F = 100((max(Q)+1)/10) =  90%
            %if N=9, min F = 100((min(Q)+1)/10) =  10%
            f = exceedance_probability(Q);
            assertEqual(testCase,min(f.F), 10);
            assertEqual(testCase,max(f.F), 90);
        end

        function test_polynomial_fit(testCase)
            % Calculate a first order polynomial on an x=y line
            poly = polynomial_fit(0:1:7, 0:1:7, 1);
            % intercept should be 0
            assertEqual(testCase,poly.coef(2), 0.0,'AbsTol',0.01);
            % slope should be 1
            assertEqual(testCase,poly.coef(1), 1.0,'AbsTol',0.01);
            % r-squared should be perfect
            assertEqual(testCase,poly.fit, 1.0,'AbsTol',0.01);
        end

        function test_discharge_to_velocity(testCase)
            % Create arbitrary discharge between 0 and 8(N=9)
            Q = struct('Discharge',[0;1;2;3;4;5;6;7;8],'time',[0 1 2 3 4 5 6 7 8]);
            % Calculate a first order polynomial on an DV_Curve x=y line 10 times greater than the Q values
            poly = polynomial_fit(0:1:8, 10*(0:1:8), 1);
            % Becuase the polynomial line fits perfect we should expect the V to equal 10*Q
            V = discharge_to_velocity(Q,poly.coef);
            assertEqual(testCase,sum(10*Q.Discharge - V.V), 0.00,'AbsTol',0.01);
        end

        function test_velocity_to_power(testCase)
            % Calculate a first order polynomial on an DV_Curve x=y line 10 times greater than the Q values
            poly = polynomial_fit(0:1:8, 10*(0:1:8),1);
%           Becuase the polynomial line fits perfect we should expect the V to equal 10*Q
            Q = struct('Discharge',[0;1;2;3;4;5;6;7;8],'time',[0 1 2 3 4 5 6 7 8]);
            V = discharge_to_velocity(Q, poly.coef);
            VV = V.V;
            % Calculate a first order polynomial on an VP_Curve x=y line 10 times greater than the V values
            poly2 = polynomial_fit(0:1:8, 10*(0:1:8),1);
            % Set cut in/out to exclude 1 bin on either end of V range
            cut_in  = VV(1);
            cut_out = VV(end);
            % Power should be 10x greater and exclude the ends of V
            P = velocity_to_power(V, poly.coef, cut_in, cut_out);
            %Cut in power zero
            assertEqual(testCase,P.P(1), 0.00,'AbsTol',0.01);
            %Cut out power zero
            assertEqual(testCase,P.P(end), 800,'AbsTol',0.01);
            % Middle 10x greater than velocity
            assertEqual(testCase,sum((P.P - 10*V.V)), 0.00,'AbsTol',0.1);
        end

        function test_energy_produced(testCase)
            seednum = 123;
            rng(seednum);
            % If power is always X then energy produced with be x*seconds 
            X=1;
            seconds=1;
            P = struct('P',[X;X;X;X;X;X;X;X;X;X],'time',[0 1 2 3 4 5 6 7 8 9]);
            EP = energy_produced(P, seconds);
            assertEqual(testCase,EP, X*seconds,'AbsTol',0.1);
            % for a normal distribution of Power EP = mean *seconds
            mu=5;
            sigma=1;
            
            function normrnd = normrnd(mu, sigma)
                normrnd = randn * sigma + mu;
            end

            power_dist = struct('P', ...
                [normrnd(mu,sigma); ...
                normrnd(mu,sigma); ...
                normrnd(mu,sigma); ...
                normrnd(mu,sigma); ...
                normrnd(mu,sigma); ...
                normrnd(mu,sigma); ...
                normrnd(mu,sigma); ...
                normrnd(mu,sigma); ...
                normrnd(mu,sigma); ...
                normrnd(mu,sigma)], ...
                'time',[0 1 2 3 4 5 6 7 8 9]);
            EP2 = energy_produced(power_dist, seconds);
            assertEqual(testCase,EP2, mu*seconds,'AbsTol',0.1);
        end

        function test_plot_flow_duration_curve(testCase)
            filename = 'river_plot_flow_duration_curve.png';
            if isfile(filename)
                delete(filename);
            end
                
            Q = struct('Discharge',[0;1;2;3;4;5;6;7;8],'time',[0 1 2 3 4 5 6 7 8]);
            f = exceedance_probability(Q);
            
            plot_flow_duration_curve(Q.Discharge, f.F,"savepath",filename);
            
            assertTrue(testCase,isfile(filename));
            delete(filename);
        end
            
        function test_plot_power_duration_curve(testCase)
            filename = 'river_plot_power_duration_curve.png';
            if isfile(filename)
                delete(filename);
            end
            
            Q = struct('Discharge',[0;1;2;3;4;5;6;7;8],'time',[0 1 2 3 4 5 6 7 8]);
            f = exceedance_probability(Q);
            poly = polynomial_fit(0:1:8, 10*(0:1:8),1);
            V = discharge_to_velocity(Q, poly.coef);
            VV = V.V;
            % Calculate a first order polynomial on an VP_Curve x=y line 10 times greater than the V values
            poly2 = polynomial_fit(0:1:8, 10*(0:1:8),1);
            % Set cut in/out to exclude 1 bin on either end of V range
            cut_in  = VV(1);
            cut_out = VV(end);
            % Power should be 10x greater and exclude the ends of V
            P = velocity_to_power(V, poly.coef, cut_in, cut_out);
            
            plot_power_duration_curve(P.P, f.F,"savepath",filename);
            
            assertTrue(testCase,isfile(filename));
            delete(filename);
        end
            
        function test_plot_velocity_duration_curve(testCase)
            filename = 'river_plot_velocity_duration_curve.png';
            if isfile(filename)
                delete(filename);
            end
            
            Q = struct('Discharge',[0;1;2;3;4;5;6;7;8],'time',[0 1 2 3 4 5 6 7 8]);
            poly = polynomial_fit(0:1:8, 10*(0:1:8),1);
            V = discharge_to_velocity(Q, poly.coef);
            f = exceedance_probability(Q);
            
            plot_velocity_duration_curve(V.V, f.F,"savepath",filename);
            
            assertTrue(testCase,isfile(filename));
            delete(filename);
        end
        
        function test_plot_discharge_timeseries(testCase)
            filename = 'river_plot_discharge_timeseries.png';
            if isfile(filename)
                delete(filename);
            end
            
            Q = struct('Discharge',[0;1;2;3;4;5;6;7;8],'time',[0 1 2 3 4 5 6 7 8]);

            plot_discharge_timeseries(Q,"savepath",filename);
            
            assertTrue(testCase,isfile(filename));
            delete(filename);
        end
            
        function test_plot_discharge_vs_velocity(testCase)
            filename = 'river_plot_discharge_vs_velocity.png';
            if isfile(filename)
                delete(filename);
            end
            
            Q = struct('Discharge',[0;1;2;3;4;5;6;7;8],'time',[0 1 2 3 4 5 6 7 8]);
            poly = polynomial_fit(0:1:8, 10*(0:1:8),1);
            V = discharge_to_velocity(Q, poly.coef);
            
            plot_discharge_vs_velocity(Q.Discharge,V.V,"savepath",filename);
            
            assertTrue(testCase,isfile(filename));
            delete(filename);
        end
        
        function test_plot_velocity_vs_power(testCase)
            filename = 'river_plot_velocity_vs_power.png';
            if isfile(filename)
                delete(filename);
            end
            
            Q = struct('Discharge',[0;1;2;3;4;5;6;7;8],'time',[0 1 2 3 4 5 6 7 8]);
            f = exceedance_probability(Q);
            poly = polynomial_fit(0:1:8, 10*(0:1:8),1);
            V = discharge_to_velocity(Q, poly.coef);
            VV = V.V;
            % Calculate a first order polynomial on an VP_Curve x=y line 10 times greater than the V values
            poly2 = polynomial_fit(0:1:8, 10*(0:1:8),1);
            % Set cut in/out to exclude 1 bin on either end of V range
            cut_in  = VV(1);
            cut_out = VV(end);
            % Power should be 10x greater and exclude the ends of V
            P = velocity_to_power(V, poly.coef, cut_in, cut_out);
            
            plot_velocity_vs_power(V.V, P.P,"polynomial_coeff",poly.coef,"savepath",filename);

            assertTrue(testCase,isfile(filename));
            delete(filename);
            end
        
    end
end

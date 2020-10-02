classdef River_TestDevice < matlab.unittest.TestCase

    methods (Test) 

        function test_circular(testCase)
            Obj.diameter = 1;

            [eq, ca] = circular(Obj.diameter);
            assertEqual(testCase,eq, Obj.diameter);
            assertEqual(testCase,ca, 4*pi*Obj.diameter^2.);
        end

        function test_ducted(testCase)
            Obj.diameter = 1;

            [eq, ca] = ducted(Obj.diameter); 
            assertEqual(testCase,eq, Obj.diameter);
            assertEqual(testCase,ca, 4*pi*Obj.diameter^2.);
        end

        function test_rectangular(testCase)
            Obj.height = 2;
            Obj.width = 3;

            [eq, ca] = rectangular(Obj.height, Obj.width);
            assertEqual(testCase,eq, 2.76, 'AbsTol',0.01)
            assertEqual(testCase,ca, Obj.height*Obj.width,'AbsTol', 0.01)
        end

        function test_multiple_circular(testCase)
            Obj.diameters = [1,2,3,4];

            [eq, ca] = multiple_circular(Obj.diameters);
            assertEqual(testCase,eq, 5.48,'AbsTol',0.01);
            assertEqual(testCase,ca, 23.56,'AbsTol', 0.01);
        end
        
        function test_tip_speed_ratio(testCase)
            rotor_speed = [15,16,17,18]; % create array of rotor speeds
            rotor_diameter = 77; % diameter of rotor for GE 1.5
            inflow_speed = [13,13,13,13]; % array of wind speeds
            TSR_answer = [4.7,5.0,5.3,5.6];
        
            TSR = tip_speed_ratio(rotor_speed/60,rotor_diameter,inflow_speed);
            assertEqual(testCase,TSR, TSR_answer,'AbsTol',0.05);
            
        end
        
        function test_power_coefficient(testCase)
            inflow_speed = [4,6,8,10,12,14,16,18,20];
            power_out = [59,304,742,1200,1400,1482,1497,1497,1511];
            capture_area = 4656.63;
            rho = 1.225;
            Cp_answer = [0.320,0.493,0.508,0.421,0.284,0.189,0.128,0.090,0.066];
        
            Cp = power_coefficient(power_out*1000,inflow_speed,capture_area,rho);
            assertEqual(testCase,Cp, Cp_answer,'AbsTol',0.01);
            
        end
    end
end
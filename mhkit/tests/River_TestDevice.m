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
    end
end
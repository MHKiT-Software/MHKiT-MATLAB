classdef Mooring_TestMooring < matlab.unittest.TestCase
    properties
        dfile
        ifile
    end
    methods (TestMethodSetup)
        function setupTestData(testCase)
            testCase.dfile = "../../examples/data/mooring/line1_test.out";
            testCase.ifile = "../../examples/data/mooring/TestInput.MD.dat";
        end

    end

    methods (Test)
        % Test methods

        function test_read_moordyn(testCase)
            [data, input] = read_moordyn(testCase.dfile, testCase.ifile);

            assertClass(testCase,data,'table')
            assertClass(testCase,input,'struct')
        end

        function test_lay_length(testCase)
            data = read_moordyn(testCase.dfile);

            ll = lay_length(data, -56, 0.25);
            ll_avg = mean(ll);
            answer = 45.0;

            assertEqual(testCase, ll_avg, answer, 'AbsTol', 0.01)
        end

        function test_plot_animation(testCase)
            data = read_moordyn(testCase.dfile);
            % test 2d
            plot_mooring_animate(data(data.Time<10, :), '2D', 0.001, ...
                'xlabel','X-axis','ylabel','depth [m]','title','2D example');
            saveas(gcf,'test_mooring_2d.png')
            % test 3d
            plot_mooring_animate(data(data.Time<10, :), '3D', 0.001, ...
                'xlabel','X-axis','ylabel','depth [m]','title','3D example');
            saveas(gcf,'test_mooring_3d.png')
            
            assertTrue(testCase,isfile('test_mooring_2d.png'))
            delete('test_mooring_2d.png')
            assertTrue(testCase,isfile('test_mooring_3d.png'))
            delete('test_mooring_3d.png')
        end
    end

end
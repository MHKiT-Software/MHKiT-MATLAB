classdef Utils_TestDefineBins < matlab.unittest.TestCase

    methods (Test)

        function test_mhkit_define_bins_1d(testCase)
            % One shared implementation for all three spec forms
            u = mhkit_define_bins_1d(struct('start', 0, 'stop', 2, 'width', 0.5));
            e = mhkit_define_bins_1d(struct('edges', 0:0.5:2));
            c = mhkit_define_bins_1d(struct('centers', 0.25:0.5:1.75));
            assertEqual(testCase, u.edges, (0:0.5:2)', 'AbsTol', 1e-12);
            assertEqual(testCase, u.centers, (0.25:0.5:1.75)', 'AbsTol', 1e-12);
            assertEqual(testCase, e, u, 'AbsTol', 1e-12);
            assertEqual(testCase, c, u, 'AbsTol', 1e-12);
            err_id = 'MHKiT:mhkit_define_bins_1d:InvalidInput';
            testCase.verifyError(@() mhkit_define_bins_1d(struct('start', 0, 'stop', 1)), err_id);
            testCase.verifyError(@() mhkit_define_bins_1d(struct('start', 0, 'stop', 3.7, 'width', 0.5)), err_id);
            testCase.verifyError(@() mhkit_define_bins_1d(struct('edges', [0 1], 'centers', 0.5)), err_id);
        end

        function test_mhkit_define_bins_2d(testCase)
            % Uniform grid from start, stop, width
            bins = mhkit_define_bins_2d(struct('start', 0, 'stop', 2, 'width', 0.5), ...
                                        struct('start', 0, 'stop', 3, 'width', 1));
            assertEqual(testCase, bins.x.centers, [0.25; 0.75; 1.25; 1.75], 'AbsTol', 1e-12);
            assertEqual(testCase, bins.x.edges, (0:0.5:2)', 'AbsTol', 1e-12);
            assertEqual(testCase, bins.y.centers, [0.5; 1.5; 2.5], 'AbsTol', 1e-12);
            assertEqual(testCase, bins.y.edges, (0:3)', 'AbsTol', 1e-12);
            % Floating point range that is a whole number of widths stays exact
            b = mhkit_define_bins_2d(struct('start', 0, 'stop', 0.1 * 3, 'width', 0.1), struct('start', 0, 'stop', 1, 'width', 1));
            assertEqual(testCase, numel(b.x.centers), 3);
            % Edges in: centers are midpoints. Centers in: edges extend half a
            % spacing beyond the end centers. Non-uniform spacing is allowed.
            c = mhkit_define_bins_2d(struct('edges', [0 0.5 1 2]), struct('centers', [0.5 1.5 2.5]));
            assertEqual(testCase, c.x.edges, [0; 0.5; 1; 2]);
            assertEqual(testCase, c.x.centers, [0.25; 0.75; 1.5]);
            assertEqual(testCase, c.y.centers, [0.5; 1.5; 2.5]);
            assertEqual(testCase, c.y.edges, [0; 1; 2; 3]);
            % Same grid as the uniform spec when the spacing is uniform
            u = mhkit_define_bins_2d(struct('start', 0, 'stop', 3, 'width', 1), struct('start', 0, 'stop', 3, 'width', 1));
            assertEqual(testCase, c.y, u.y, 'AbsTol', 1e-12);
            % Every error carries this function's name, even those raised
            % inside the shared mhkit_define_bins_1d
            err_id = 'MHKiT:mhkit_define_bins_2d:InvalidInput';
            ok = struct('start', 0, 'stop', 1, 'width', 1);
            testCase.verifyError(@() mhkit_define_bins_2d(struct('start', 2, 'stop', 1, 'width', 0.5), ok), err_id);
            testCase.verifyError(@() mhkit_define_bins_2d(struct('start', 0, 'stop', 1, 'width', -1), ok), err_id);
            testCase.verifyError(@() mhkit_define_bins_2d(struct('start', 0, 'stop', 2.2, 'width', 1), ok), err_id);
            testCase.verifyError(@() mhkit_define_bins_2d(struct('start', 0, 'end', 1, 'width', 1), ok), err_id);
            testCase.verifyError(@() mhkit_define_bins_2d(struct('start', 0, 'stop', 1, 'width', 1, 'extra', 1), ok), err_id);
            testCase.verifyError(@() mhkit_define_bins_2d(struct('edges', [0 1], 'centers', 0.5), ok), err_id);
            testCase.verifyError(@() mhkit_define_bins_2d(struct('width', 1), ok), err_id);
            testCase.verifyError(@() mhkit_define_bins_2d(struct('edges', [1 1 2]), ok), err_id);
            testCase.verifyError(@() mhkit_define_bins_2d(struct('centers', 1), ok), err_id);
        end

    end

end

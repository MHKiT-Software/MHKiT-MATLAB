classdef Examples_TestMatlabOnly < matlab.unittest.TestCase
    % Executes the example live scripts (examples/*.mlx) that use MATLAB only code
    %
    % Examples not included because they call Python bridge code:
    %   *_contour*, extreme_response_full_sea_state, short_term_extremes,
    %   loads, power, qc, tidal, SWAN, delft3d, wecsim_example
    % Examples not included because they make live network requests:
    %   WPTO_hindcast

    properties (TestParameter)
        exampleName = { ...
            'acoustics_example', ...
            'adcp_example', ...
            'adv_example', ...
            'cdip_example', ...
            'mooring_example', ...
            'power_example', ...
            'river_example', ...
            'strain_measurement_example', ...
            'upcrossing_example', ...
            'wave_example', ...
            'wecsim_power_performance_example', ...
        };
    end

    properties (Constant)
        repoRoot = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
        examplesFolder = fullfile(Examples_TestMatlabOnly.repoRoot, 'examples');
        % Examples that open data files relative to the repository root
        % (e.g. './examples/data/...'). All others use paths relative to the
        % examples folder (e.g. './data/...').
        runFromRepoRoot = {'mooring_example'};
    end

    methods (TestMethodTeardown)
        function closeFigures(~)
            close all force;
        end
    end

    methods (Test)
        function test_example_runs(testCase, exampleName)
            import matlab.unittest.fixtures.CurrentFolderFixture
            import matlab.unittest.fixtures.PathFixture

            file = fullfile(testCase.examplesFolder, [exampleName '.mlx']);
            testCase.assertTrue(isfile(file), ['Example not found: ' file]);

            if ismember(exampleName, testCase.runFromRepoRoot)
                workingFolder = testCase.repoRoot;
            else
                workingFolder = testCase.examplesFolder;
            end
            % Fixtures restore the path and folder when the test finishes.
            testCase.applyFixture(PathFixture(testCase.examplesFolder));
            testCase.applyFixture(CurrentFolderFixture(workingFolder));

            runInIsolatedWorkspace(exampleName);
        end
    end
end

function runInIsolatedWorkspace(exampleName)
    % Runs the live script by name (not run(path), which would change the
    % current folder) in its own workspace, so 'clear' inside an example
    % cannot touch the test framework. Output is captured to keep logs short.
    evalc(exampleName);
end

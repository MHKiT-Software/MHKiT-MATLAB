function results = runMatlabOnlyExamples()
% runMatlabOnlyExamples Execute the MHKiT example live scripts built from native MATLAB code
%

    import matlab.unittest.TestSuite
    import matlab.unittest.TestRunner

    testsFolder = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(fileparts(testsFolder));
    addpath(genpath(fullfile(repoRoot, 'mhkit')));

    suite = TestSuite.fromFile(fullfile(testsFolder, 'examples', 'Examples_TestMatlabOnly.m'));
    fprintf('Running %d MATLAB-only examples\n', numel(suite));

    runner = TestRunner.withTextOutput;
    results = runner.run(suite);
    disp(table(results));
end

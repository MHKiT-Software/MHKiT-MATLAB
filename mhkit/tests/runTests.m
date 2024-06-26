function results = runTests()

    import matlab.unittest.TestRunner;
    import matlab.unittest.TestSuite;
    import matlab.unittest.plugins.TestReportPlugin;
    import matlab.unittest.plugins.CodeCoveragePlugin;
    import matlab.unittest.plugins.codecoverage.CoverageReport;


    % Define test suite
    testFileName = mfilename('fullpath');
    testsFolder = fileparts(testFileName);
    suite = TestSuite.fromFolder(testsFolder);

    % Build the runner
    runner = TestRunner.withTextOutput;

    % Add HTML plugin
    htmlFolder = fullfile(testsFolder, 'test_results');
    plugin = TestReportPlugin.producingHTML(htmlFolder);
    runner.addPlugin(plugin);

    % Add PDF
    pdfFile = fullfile(htmlFolder, 'test_results.pdf');
    plugin = TestReportPlugin.producingPDF(pdfFile);
    runner.addPlugin(plugin);

    % Add Code Coverage Report
    sourceCodeFolder = fileparts(testsFolder);
    codeCoverageFolder = fullfile(testsFolder, 'coverage_report');
    reportFormat = CoverageReport(codeCoverageFolder);

    p = CodeCoveragePlugin.forFolder(sourceCodeFolder, Producing=reportFormat, IncludingSubfolders=true);
    runner.addPlugin(p);

    % Run the tests
    results = runner.run(suite);

end


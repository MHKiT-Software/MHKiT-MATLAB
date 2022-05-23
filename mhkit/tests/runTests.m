function results = runTests()

    import matlab.unittest.TestRunner;
    import matlab.unittest.TestSuite;
    import matlab.unittest.plugins.TestReportPlugin;

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

    % Run the tests
    results = runner.run(suite);

end
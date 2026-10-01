function results = runPythonTests(pythonExecutable, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Run the full unit test suite, including tests that call MHKiT-Python
%
% Configures MATLAB to run the Python interpreter at pythonExecutable
% OutOfProcess, so that the BLAS, HDF5, and netCDF libraries loaded by
% Python do not conflict with the versions that ship with MATLAB, checks
% that the Python modules the tests use can be imported, then runs every
% test in mhkit/tests. Run locally with
% scripts/run_python_tests_macos.sh on macOS or
% scripts/run_python_tests_windows.ps1 on Windows.
%
% On macOS the OutOfProcess Python host still resolves libexpat, libssl,
% and libcrypto from MATLAB's bin/maca64 folder. These are older than the
% versions conda builds of Python, pyexpat, and netCDF4 link against,
% which causes errors like "Symbol not found:
% _XML_SetAllocTrackerActivationThreshold" and "Unable to resolve the
% name 'py.mhkit...'". Preloading the Python environment's copies with
% DYLD_INSERT_LIBRARIES fixes this without affecting MATLAB.
%
% On Windows the same conflict appears as "DLL load failed while importing
% pyexpat" with conda builds of Python, which load a shared libexpat.dll.
% MATLAB's copy is already loaded, so changing PATH cannot fix it. Use
% Python from python.org, which builds expat into pyexpat, with packages
% from PyPI. The Python environment's folders are still placed first on
% PATH so its other DLLs are found.
%
% Parameters
% ------------
% pythonExecutable : string
%   Path to the Python executable. The environment must have mhkit and
%   mhkit_python_utils installed.
% ImportCheckOnly : logical (optional)
%   Name-value argument. Only check that the Python modules the tests use
%   can be imported, without running the tests. Default false.
%
% Returns
% ---------
% results : matlab.unittest.TestResult array
%   Results for every test in mhkit/tests
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    arguments
        pythonExecutable {mustBeTextScalar}
        options.ImportCheckOnly (1,1) logical = false
    end

    import matlab.unittest.TestSuite
    import matlab.unittest.TestRunner

    if ismac
        envLib = fullfile(fileparts(fileparts(char(pythonExecutable))), 'lib');
        libs = fullfile(envLib, {'libexpat.1.dylib', 'libssl.3.dylib', 'libcrypto.3.dylib'});
        missing = libs(~isfile(libs));
        if ~isempty(missing)
            warning('MHKiT:runPythonTests:MissingLibrary', ...
                ['Python environment libraries not found, MATLAB''s bundled versions ' ...
                'will be used instead and Python imports may fail: %s'], strjoin(missing, ', '));
        end
        libs = libs(isfile(libs));
        setenv('DYLD_INSERT_LIBRARIES', strjoin(libs, ':'));
        fprintf('DYLD_INSERT_LIBRARIES=%s\n', getenv('DYLD_INSERT_LIBRARIES'));
    elseif ispc
        envRoot = fileparts(char(pythonExecutable));
        envDirs = fullfile(envRoot, {'', 'Library\mingw-w64\bin', 'Library\usr\bin', 'Library\bin', 'Scripts'});
        envDirs = envDirs(isfolder(envDirs));
        setenv('PATH', strjoin([envDirs, {getenv('PATH')}], pathsep));
        fprintf('PATH=%s\n', getenv('PATH'));
    end

    pe = pyenv(Version=pythonExecutable, ExecutionMode="OutOfProcess");
    disp(pe);

    check_python_imports();
    if options.ImportCheckOnly
        results = matlab.unittest.TestResult.empty;
        return
    end

    testsFolder = fileparts(mfilename('fullpath'));
    repoRoot = fileparts(fileparts(testsFolder));
    addpath(genpath(fullfile(repoRoot, 'mhkit')));

    suite = TestSuite.fromFolder(testsFolder);
    runner = TestRunner.withTextOutput;
    results = runner.run(suite);
end

function check_python_imports()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check that the Python modules the tests depend on can be imported
%
% Fails with one clear error if any module is missing or cannot load its
% native libraries. Without this check these failures appear as "Unable
% to resolve the name 'py.mhkit...'" in every test that calls
% MHKiT-Python.
%
% Parameters
% ------------
% None
%
% Returns
% ---------
% None
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    modules = {'pyexpat', 'ssl', 'numpy', 'scipy', 'pandas', 'xarray', ...
        'netCDF4', 'h5py', 'pecos', 'mhkit', 'mhkit.loads', 'mhkit.qc', ...
        'mhkit.river', 'mhkit.wave', 'mhkit_python_utils'};

    failures = {};
    for i = 1:numel(modules)
        try
            module = py.importlib.import_module(modules{i});
            % Without MHKiT-Python installed, the MATLAB mhkit folder in the
            % current folder imports as an empty namespace package
            if isa(py.getattr(module, '__file__', py.None), 'py.NoneType')
                error('MHKiT:runPythonTests:NamespacePackage', ...
                    'imported an empty namespace package from %s, the module is not installed', ...
                    char(py.str(py.getattr(module, '__path__', py.None))));
            end
            fprintf('Python import OK:     %s\n', modules{i});
        catch err
            fprintf('Python import FAILED: %s\n', modules{i});
            failures{end+1} = sprintf('  %s: %s', modules{i}, err.message); %#ok<AGROW>
        end
    end

    % Versions of the native libraries Python loaded, MATLAB's bundled
    % versions appearing here indicates a library conflict
    try
        fprintf('Python version:  %s\n', char(py.sys.version));
        fprintf('expat version:   %s\n', char(py.pyexpat.EXPAT_VERSION));
        fprintf('OpenSSL version: %s\n', char(py.ssl.OPENSSL_VERSION));
        fprintf('MHKiT-Python:    %s\n', char(py.getattr(py.importlib.import_module('mhkit'), '__version__')));
    catch
    end

    if ~isempty(failures)
        error('MHKiT:runPythonTests:ImportFailed', ...
            ['%d Python module(s) failed to import, so tests that call ' ...
            'MHKiT-Python cannot run:\n%s\n\nIf the errors mention symbols or ' ...
            'DLLs (pyexpat, libssl, libcurl), MATLAB''s bundled libraries are ' ...
            'being loaded instead of the Python environment''s. See the ' ...
            'runPythonTests help for details.'], numel(failures), strjoin(failures, newline));
    end
end

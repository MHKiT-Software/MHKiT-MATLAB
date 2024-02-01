function results = runCodeCompatibilityReport()
    % Find source code folder
    this_filename = mfilename('fullpath');
    this_folder = fileparts(this_filename);
    sourceCodeFolder = fileparts(this_folder);

    % Generate MATLAB code compatibility report
    r = analyzeCodeCompatibility(sourceCodeFolder, 'IncludeSubfolders', true);

    % Print code compatability report
    display("MHKiT MATLAB Code Compatibility Report");
    display(r.Recommendations);

    % Search for error string in the `Severity` column. If found exit with status code 1 to indicate to the test runner that ther
    error_string = 'Error';
    severity_array = r.Recommendations.Severity;

    % Convert categorical array to cell array of strings
    severity_cell_array = cellstr(severity_array);

    % Sanitize error string
    sanitized_error_string = lower(strtrim(error_string));

    % Sanitize each element in the array using cellfun
    sanitized_severity_array = cellfun(@(x) lower(strtrim(x)), severity_cell_array, 'UniformOutput', false);

    if any(strcmp(sanitized_severity_array, sanitized_error_string))
        disp(['Code Compatability Failure: Found "', error_string, '". in the MHKiT-MATLAB code compatability report. See above output for details.']);
        exit(1);
    else
        disp('Success: MHKiT-MATLAB Code Compatability Analysis did not find errors!');
        exit(0);
    end


    results = runCodeCompatibilityReport();

end

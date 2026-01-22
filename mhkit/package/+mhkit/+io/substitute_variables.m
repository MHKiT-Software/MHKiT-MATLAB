function result = substitute_variables(text, spec)
    % Substitute variables in text using spec values (shared utility)
    %
    % Parameters:
    %   text (string): Text with variables to substitute
    %   spec (struct): Specification structure
    %
    % Returns:
    %   result (string): Text with variables substituted

    result = text;

    % Standard substitutions
    if contains(result, '<conda_env>')
        result = strrep(result, '<conda_env>', spec.conda.environment_name);
    end

    if contains(result, '<python_version>')
        result = strrep(result, '<python_version>', spec.python.install_version);
    end

    if contains(result, '<mhkit_python_version>')
        result = strrep(result, '<mhkit_python_version>', spec.mhkit_python.version);
    end

    % Platform-specific substitutions
    if contains(result, '<conda_lib_path>')
        conda_env = spec.conda.environment_name;
        python_cmd = mhkit.sys.python_cmd();
        [status, conda_info] = mhkit.sys(sprintf('conda run -n %s %s -c "import sys; import os; print(os.path.join(os.path.dirname(sys.executable), ''lib''))"', conda_env, python_cmd));
        if status == 0
            conda_lib_path = strtrim(conda_info);
            result = strrep(result, '<conda_lib_path>', conda_lib_path);
        end
    end
end
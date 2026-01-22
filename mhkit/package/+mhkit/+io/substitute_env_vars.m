function result = substitute_env_vars(template)
    % Substitute environment variables in template
    %
    % Parameters:
    %   template (string): Template with {ENV_VAR} placeholders
    %
    % Returns:
    %   result (string): Template with environment variables substituted

    result = template;

    % Replace {USERPROFILE} with actual value
    if contains(result, '{USERPROFILE}')
        result = strrep(result, '{USERPROFILE}', getenv('USERPROFILE'));
    end

    % Replace {HOME} with actual value
    if contains(result, '{HOME}')
        result = strrep(result, '{HOME}', getenv('HOME'));
    end

    % Convert forward slashes to platform-appropriate separators
    result = strrep(result, '/', filesep);
end
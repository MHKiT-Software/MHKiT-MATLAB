function [status, cmdout] = sys(command)
    % SYS_SHELL Cross-platform initialized shell execution
    % Executes commands in a shell with the user's environment properly initialized
    %
    % Inputs:
    %   command: String containing the command to execute
    %
    % Outputs:
    %   status: Exit status (0 for success)
    %   cmdout: Command output as string

    % Clean command by escaping quotes appropriately for the target shell
    if ispc
        % Windows: double quotes for cmd.exe
        cleaned_command = strrep(command, '"', '""');
        % Windows: Use cmd.exe with user's environment
        full_cmd = sprintf('cmd.exe /c "%s"', cleaned_command);
    else
        % Mac/Linux: escape double quotes with backslash
        cleaned_command = strrep(command, '"', '\"');
        % Mac/Linux: Use login shell to get full environment
        full_cmd = sprintf('bash -l -c "%s"', cleaned_command);
    end

    disp("Executing...");
    disp(full_cmd);

    [status, cmdout] = system(full_cmd);
end

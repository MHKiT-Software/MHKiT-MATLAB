function [status, cmdout] = sys(command)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Cross-platform initialized shell execution
%   Executes commands in a shell with the user's environment properly initialized
%
% Parameters
% ------------
%     command: string
%         Command to execute in the shell
%
% Returns
% ---------
%     status: double
%         Exit status (0 for success, non-zero for error)
%
%     cmdout: string  
%         Command output as string
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    [status, cmdout] = system(full_cmd);
end

function cmd = python_cmd()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Get platform-appropriate Python command
%
% Parameters
% ------------
%     None
%
% Returns
% ---------
%     cmd: string
%         Python command ('python.exe' on Windows, 'python' elsewhere)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ispc
        cmd = 'python.exe';
    else
        cmd = 'python';
    end
end
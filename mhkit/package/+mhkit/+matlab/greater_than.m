function result = greater_than(maximum_version)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check if current MATLAB version is greater than specified version
%
% Parameters
% ------------
%     maximum_version: string
%         Maximum tested MATLAB version (e.g., 'R2024b')
%
% Returns
% ---------
%     result: logical
%         true if current MATLAB version is greater than maximum_version
%         false otherwise
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check if current MATLAB version is newer than the specified version
    % isMATLABReleaseOlderThan returns true if current < specified
    % We want true if current > specified
    current_release = version('-release');
    result = ~isMATLABReleaseOlderThan(maximum_version) && ~strcmp(current_release, maximum_version);
end

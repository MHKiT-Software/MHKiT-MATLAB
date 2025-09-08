function result = less_than(minimum_version)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Check if current MATLAB version is less than specified version
%
% Parameters
% ------------
%     minimum_version: string
%         Minimum required MATLAB version (e.g., 'R2022b')
%
% Returns
% ---------
%     result: logical
%         true if current MATLAB version is less than minimum_version
%         false otherwise
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    result = isMATLABReleaseOlderThan(minimum_version);
end

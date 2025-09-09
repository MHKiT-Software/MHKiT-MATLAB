function platform = get_platform()
    % Get platform string for hook execution
    %
    % Returns:
    %   platform (string): Platform identifier ('windows', 'linux', 'mac')
    
    if ispc
        platform = 'windows';
    elseif ismac
        platform = 'mac';
    elseif isunix
        platform = 'linux';
    else
        platform = 'unknown';
    end
end
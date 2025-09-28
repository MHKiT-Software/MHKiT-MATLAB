function choice = get_user_consent(message, code_preview)
    % Get user consent with standardized y/n/exit options
    %
    % Parameters:
    %   message (string): Message to display to user
    %   code_preview (cell array, optional): Code lines to display
    %
    % Returns:
    %   choice (string): User choice ('y', 'n', or 'exit')

    if nargin < 2
        code_preview = {};
    end

    fprintf('\n%s\n', message);

    % Display code preview if provided
    if ~isempty(code_preview)
        for i = 1:length(code_preview)
            fprintf('  %s\n', code_preview{i});
        end
        fprintf('\n');
    end

    % Get user input with validation
    while true
        choice = input('', 's');
        choice = lower(strtrim(choice));

        if ismember(choice, {'y', 'yes'})
            choice = 'y';
            break;
        elseif ismember(choice, {'n', 'no'})
            choice = 'n';
            break;
        elseif ismember(choice, {'exit', 'quit', 'e', 'q'})
            choice = 'exit';
            break;
        else
            fprintf('Please enter y, n, or exit: ');
        end
    end
end
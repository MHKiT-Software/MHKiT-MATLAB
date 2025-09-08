function logger = get_logger()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Returns a simple logger struct with basic logging functions
%
% Parameters
% ------------
%     None
%
% Returns
% ---------
%     logger: struct
%         Logger structure with debug, info, warning, and error functions
%         Each function accepts a format string and optional arguments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    logger = struct();
    logger.debug = @(msg, varargin) fprintf(['MHKiT Debug: ' msg '\n'], varargin{:});
    logger.info = @(msg, varargin) fprintf('%s\n', sprintf(msg, varargin{:}));
    logger.warning = @(msg, varargin) fprintf('MHKiT Warning: %s\n', sprintf(msg, varargin{:}));
    logger.error = @(msg, varargin) fprintf('MHKiT Error: %s\n', sprintf(msg, varargin{:}));
end

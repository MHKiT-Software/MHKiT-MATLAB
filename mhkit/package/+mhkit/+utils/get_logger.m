function logger = get_logger()
    % GET_LOGGER Returns a simple logger struct with basic logging functions
    logger = struct();
    logger.debug = @(msg, varargin) fprintf(['MHKiT Debug: ' msg '\n'], varargin{:});
    logger.info = @(msg, varargin) fprintf('MHKiT: %s\n', sprintf(msg, varargin{:}));
    logger.warning = @(msg, varargin) fprintf('MHKiT Warning: %s\n', sprintf(msg, varargin{:}));
    logger.error = @(msg, varargin) fprintf('MHKiT Error: %s\n', sprintf(msg, varargin{:}));
end

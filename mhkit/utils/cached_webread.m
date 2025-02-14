function data = cached_webread(url, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Perform webread with caching to avoid redundant downloads.
%
%     Parameters
%     ----------
%         url : char
%             The URL from which to retrieve the data.
%
%         options : weboptions
%             (Optional) Additional options for the webread function.
%
%      Returns
%      -------
%         data : <Url Dependent (json, XML, plain text, etc)>
%             The data retrieved from the specified URL.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Configuration
    maxRetries = 5;
    retryPause = 5;

    % Setup cache directory
    cacheDir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'mhkit_webread_cache');
    if ~exist(cacheDir, 'dir')
        mkdir(cacheDir);
    end

    % Create URL hash for cache
    cacheFilename = fullfile(cacheDir, [md5hash(url) '.mat']);

    % Check cache
    if exist(cacheFilename, 'file')
        cachedFileinfo = dir(cacheFilename);
        cachedFileAge = datetime('now') - datetime(cachedFileinfo.datenum, 'ConvertFrom', 'datenum');

        if days(cachedFileAge) <= 1
            fprintf('Reading from cache: %s\n', truncateUrl(url));
            load(cacheFilename, 'data');
            return;
        end
    end

    % Check file size before downloading
    try
        headOptions = weboptions('RequestMethod', 'head');
        response = webwrite(url, '', headOptions);

        % Get content length from headers if available
        if isfield(response, 'Content-Length')
            content_length = str2double(response.('Content-Length'));
            if ~isnan(content_length)
                fprintf('Downloading %.1f MB: %s\n', content_length/1024/1024, truncateUrl(url));
            else
                fprintf('Downloading: %s\n', truncateUrl(url));
            end
        else
            fprintf('Downloading: %s\n', truncateUrl(url));
        end
    catch
        % If HEAD request fails, continue without size information
        fprintf('Downloading: %s\n', truncateUrl(url));
    end

    % Perform webread with progress tracking
    for attempt = 1:maxRetries
        try
            data = webread(url, options);
            fprintf('Successfully downloaded data\n');
            break;
        catch exception
            if attempt < maxRetries
                fprintf('Download failed, attempt %d of %d. Retrying...\n', attempt, maxRetries);
                pause(retryPause);
            else
                error('Failed to download after %d attempts', maxRetries);
            end
        end
    end

    % Cache the successful response
    save(cacheFilename, 'data');
end

function shortened_url = truncateUrl(url)
    max_length = 80;
    if length(url) > max_length
        shortened_url = [url(1:max_length/2-3) '...' url(end-max_length/2+4:end)];
    else
        shortened_url = url;
    end
end

function hash = md5hash(str)
    md = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(str));
    hashBytes = typecast(md.digest, 'uint8');
    hash = sprintf('%02x', hashBytes);
end

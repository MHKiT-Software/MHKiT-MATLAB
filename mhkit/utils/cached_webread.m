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

    % Define the fixed cache directory at the project root
    cacheDir = fullfile(fileparts(mfilename('fullpath')), '..', '..' , 'mhkit_webread_cache');

    % Create the cache directory if it does not exist
    if ~exist(cacheDir, 'dir')
        mkdir(cacheDir);
    end

    % Create a hash of the URL to use as a cache filename
    cacheFilename = fullfile(cacheDir, [md5hash(url) '.mat']);

    % Check if the file is in the cache
    if exist(cacheFilename, 'file')
        % Check the age of the cache
        cachedFileinfo = dir(cacheFilename);
        cachedFileAge = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss') - datetime(cachedFileinfo.datenum, 'ConvertFrom', 'datenum');

        if days(cachedFileAge) <= 1
            disp(['Data found in cache. Using cached version of ', url]);
            load(cacheFilename, 'data');
            return;
        else
            disp(['Cache is older than 1 day. Redownloading...']);
        end
    else
        % disp(['Cache not found or forceAddDate is true. Redownloading...']);
    end

    % Perform webread
    data = webread(url, options);

    % Save raw data to cache
    save(cacheFilename, 'data');
end

function hash = md5hash(str)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Compute MD5 hash for a given string.
%
%     Parameters
%     ----------
%         str : char
%             The input string for which the MD5 hash will be computed.
%
%      Returns
%      -------
%         hash : char
%             The MD5 hash value represented as a hexadecimal string.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    md = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(str));
    hashBytes = typecast(md.digest, 'uint8');
    hash = sprintf('%02x', hashBytes);
end

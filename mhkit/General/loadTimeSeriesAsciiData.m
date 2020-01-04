function [data, header] = loadTimeSeriesAsciiData(filename, varargin)
%
% [data, header] = loadTimeSeriesAsciiData(filename, delimiter, nSkipLines,filetype)
%
%
% Loads data from a ASCII text file with header lines and tranforms those
% data into a matrix. The header is saved into a separate cell array.
% Data must be in filename must be stored as either a column or row
% vectors in the datafile.
%
% Input:
%
%   filename                      name of the text data file, including
%                                 the file extension.
%   delimiter                     delimiter between data columns, default
%                                 is a space
%   nSkipLine                     skips the first nSkipLines and starts
%                                 reading data at nSkipLines+1
%
% Output:
%   data                          matrix of the data
%   header                        cell array of the header lines
%
% Dependencies
%   none
%
% Usage
% loadTextData(filename)
%   reads the CDIP data from filename and adds data to resource. Duplicate
%   times are ignored, and the output resource is sorted in ascending
%   datetime.
%
% need other usage description here %%%%
%
%
% Version 1.01, 11/25/2018Rick Driscoll, NREL

% checking to see if the file exists
if exist(filename,'file') ~= 2
    
    error(['loadTextData: ' filename ' does not exist'])
end;



switch nargin
    case 1
        % only the filename has been specified, thus spaces are assumed
        % as the delimiters
        try
            tempdata = importdata(filename);
        catch
            error(['loadTextData: unable to load data from ' filename])
        end
        data = tempdata.data;
        header = tempdata.textdata;
    case 2
        % check to see if a delimiter has been specified, otherwise, the
        % number of lines to skip has been specified
        if ischar(varargin{1})
            try
                tempdata = importdata(filename,varargin{1});
            catch
                error(['loadTextData: unable to load data from ' filename])
            end
        else 
            if isfloat(varargin{1})  & rem(varargin{1},1) == 0
                try
                    tempdata = importdata(filename,varargin{1});
                catch
                    error(['loadTextData: unable to load data from ' filename])
                end
            else
                error('loadTextData: delimiter must be a character or nSkipLines must be an integer');
            end;
        end;
        data = tempdata.data;
        header = tempdata.textdata;
    case 3
        if ~ischar(varargin{1})
            error('loadTextData: delimiter must be a character');
        end;
        if ~(isfloat(varargin{2})  & rem(varargin{2},1) == 0)
            error('loadTextData: nSkipLines must be an integer');
        end;
        try
            tempdata = importdata(filename, varargin{1}, varargin{2});
        catch
            error(['loadTextData: unable to load data from ' filename])
        end
        data = tempdata.data;
        header = tempdata.textdata;
    otherwise
        error('loadTextData: Incorrect number of input arguments')
end
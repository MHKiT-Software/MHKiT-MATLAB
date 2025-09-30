function [data, input] = read_moordyn(filepath, inputfile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads in MoorDyn OUT files such as "FAST.MD.out" and "FAST.MD.Line1.out" 
% and stores as a table. Also allows for parsing and storage of MoorDyn 
% input file as a struct.
%
%    Parameters
%    ----------
%    filepath : str
%        Path to MoorDyn OUT file
%    inputfile : str
%        Path to MoorDyn input file
%
%    Returns
%    -------
%    data : Array
%        Array containing parsed MoorDyn OUT file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    filepath {mustBeTextScalar}
    inputfile {mustBeTextScalar} = "None"
end

if nargout >= 1
    data = readtable(filepath, 'FileType','delimitedtext');
end
if nargout == 2
    if strcmp(inputfile, "None")
        error('You must provide a path to an input file!')
    end
    fid = fopen(inputfile, 'r');
    lines = {};
    while ~feof(fid)
        lines{end+1} = fgetl(fid);
    end
    fclose(fid);
    
    input = struct();
    
    % Helper function to find section start
    find_section = @(name) find(contains(lines, name), 1);
    
    % Parse Line Types
    idx = find_section('LINE TYPES');
    if ~isempty(idx)
        header = strsplit(strtrim(lines{idx+1}));
        types = {};
        i = idx+3;
        while ~isempty(lines{i}) && ~contains(lines{i}, 'POINTS')
            types{end+1} = strsplit(strtrim(lines{i}));
            i = i+1;
        end
        % Convert to table
        T = cell2table(vertcat(types{:}), 'VariableNames', header);
        % Convert numeric columns to double where possible
        for k = 2:width(T)
            T.(k) = varfun(@(x) str2double(x), T(:,k), 'OutputFormat', 'uniform');
        end
        input.LineTypes = T;
    end
    
    % Parse Points as a table
    idx = find_section('POINTS');
    if ~isempty(idx)
        header = strsplit(strtrim(lines{idx+1}));
        points = {};
        i = idx+3;
        while ~isempty(lines{i}) && ~contains(lines{i}, 'LINES')
            row = strsplit(strtrim(lines{i}));
            points(end+1, :) = row; % Collect as flat cell array
            i = i+1;
        end
        % Convert to table
        T = cell2table(points, 'VariableNames', header);
        % Convert numeric columns to double where possible (skip Attachment column)
        for k = [1,3:width(T)]
            T.(k) = str2double(T.(k));
        end
        input.Points = T;
    end
    
    % Parse Lines as a table
    idx = find_section('LINES');
    if ~isempty(idx)
        header = strsplit(strtrim(lines{idx+1}));
        linesec = {};
        i = idx+3;
        while ~isempty(lines{i}) && ~contains(lines{i}, 'SOLVER OPTIONS')
            row = strsplit(strtrim(lines{i}));
            linesec(end+1, :) = row; % Collect as flat cell array
            i = i+1;
        end
        % Convert to table
        T = cell2table(linesec, 'VariableNames', header);
        % Convert numeric columns to double where possible (skip LineType, AttachA, AttachB)
        for k = [1,4:width(T)]
            T.(k) = str2double(T.(k));
        end
        input.Lines = T;
    
    end
    % Parse Solver Options as a multiline text cell
    idx = find_section('SOLVER OPTIONS');
    if ~isempty(idx)
        solveropts = {};
        i = idx+1;
        while i <= numel(lines) && ~isempty(lines{i}) && ~contains(lines{i}, 'OUTPUTS')
            solveropts{end+1,1} = lines{i};
            i = i+1;
        end
        input.SolverOptions = solveropts;
    end
    
    % Parse Outputs
    idx = find_section('OUTPUTS');
    if ~isempty(idx)
        outputs = {};
        i = idx+1;
        while i <= numel(lines) && ~contains(lines{i}, 'END')
            outputs{end+1} = strtrim(lines{i});
            i = i+1;
        end
        input.Outputs = outputs;
    end
end




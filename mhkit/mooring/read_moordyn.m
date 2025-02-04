function data = read_moordyn(filepath)
%Reads in MoorDyn OUT files such as "FAST.MD.out" and "FAST.MD.Line1.out" 
% and stores as a table.
%    Parameters
%    ----------
%    filepath : str
%        Path to MoorDyn OUT file
%
%    Returns
%    -------
%    Array
%        Array containing parsed MoorDyn OUT file

arguments
    filepath {mustBeTextScalar}
end

data = readtable(filepath, 'FileType','delimitedtext');

end




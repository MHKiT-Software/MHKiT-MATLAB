function result = delft_3d_open_netcdf(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Opens a Delft 3D netCDF file as a python netCDF4 object
%
% Parameters
% ------------
%    filename: str
%        The name of the netCDF file to open.
%
% Returns
% ---------
%    result: py.netCDF4._netCDF4.Dataset
%        The object representing the opened netCDF dataset.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist(filename, 'file')
        ME = MException('MATLAB:delft_3d_open_netcdf:FileNotFound', "Input file '" + filename + "' does not exist.");
        throw(ME);
    end

    % Obtain a Python netCDF object. This approach leverages existing Python functionality,
    % avoiding the need to reimplement it in MATLAB, albeit requiring handling a Python object.
    result = py.mhkit_python_utils.d3d_converter.open_netcdf_file(filename);
end


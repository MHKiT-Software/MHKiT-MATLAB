import os

import netCDF4


def open_netcdf_file(filename):
    """
    Opens a Delft3D NetCDF file.

    Args:
        filename (str): The path to the NetCDF file.

    Returns:
        netCDF4.Dataset: The opened NetCDF dataset object.

    Raises:
        FileNotFoundError: If the specified file does not exist.
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"The file '{filename}' does not exist.")

    return netCDF4.Dataset(filename)

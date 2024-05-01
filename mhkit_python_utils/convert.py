import json

import xarray as xr

from datetime import date, datetime, timedelta


# https://stackoverflow.com/a/22238613
def json_serial(obj):
    """
    JSON serializer for objects not serializable by default json code.

    Parameters
    ----------
    obj : object
        Object to be serialized.

    Returns
    -------
    str
        Serialized object in ISO format or as string.

    Raises
    ------
    TypeError
        If the object type is not serializable.

    """
    if isinstance(obj, (datetime, date)):
        return obj.isoformat()
    if isinstance(obj, (datetime, timedelta)):
        return str(obj)
    raise TypeError("Type %s not serializable" % type(obj))


def xarray_to_serialized_json(ds: xr.Dataset) -> str:
    """
    Convert xarray dataset to serialized JSON.

    Parameters
    ----------
    ds : xarray.Dataset
        Dataset to convert.

    Returns
    -------
    str
        Serialized JSON representation of the dataset.

    """
    return json.dumps(ds.to_dict(), default=json_serial)


def serialized_json_to_xarray(serialized_json: str) -> xr.Dataset:
    """
    Convert serialized JSON to xarray dataset.

    Parameters
    ----------
    serialized_json : str
        Serialized JSON string representing the dataset.

    Returns
    -------
    xr.Dataset
        Xarray dataset reconstructed from the serialized JSON.

    """
    return xr.Dataset.from_dict(json.loads(serialized_json))


if __name__ == "__main__":
    from pathlib import Path
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme()

    netcdf_file = Path("../examples/data/river/d3d/turbineTest_map.nc").resolve()
    ds = xr.open_dataset(netcdf_file)
    print(ds.attrs)

    serialized_json = xarray_to_serialized_json(ds)
    ds = serialized_json_to_xarray(serialized_json)

    turkin = ds.data_vars["turkin1"]
    turkin_1d = turkin.isel(nFlowLink=list(range(1, 4)), wdim=1).plot()
    plt.show()
